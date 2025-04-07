import os
import sys
import re
import threading
import platform
from tkinter import filedialog, messagebox

import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from PIL import Image, ImageTk

from cats import run_cats

class TextRedirector:
    def __init__(self, widget):
        self.widget = widget
        self.timestamp_pattern = re.compile(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}')

    def write(self, msg):
        self.widget.configure(state='normal')

        # Colorize timestamps
        for match in self.timestamp_pattern.finditer(msg):
            start_idx, end_idx = match.span()
            # Insert everything before match as normal text
            self.widget.insert(END, msg[:start_idx], ("message",))
            # Insert matched timestamp text with a separate tag
            self.widget.insert(END, msg[start_idx:end_idx], ("timestamp",))
            msg = msg[end_idx:]

        # Colorize specific keywords
        parts = msg.split()
        for part in parts:
            if "INFO" in part:
                self.widget.insert(END, part + ' ', ("info",))
            elif "WARNING" in part:
                self.widget.insert(END, part + ' ', ("warning",))
            else:
                self.widget.insert(END, part + ' ', ("message",))

        self.widget.insert(END, '\n')
        self.widget.configure(state='disabled')
        self.widget.see(END)

    def flush(self):
        pass


def main():
    ICON_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)), "logo.png")

    def browse_file(entry):
        filename = filedialog.askopenfilename()
        entry.delete(0, END)
        entry.insert(0, filename)

    def browse_folder(entry):
        foldername = filedialog.askdirectory()
        if foldername:
            entry.delete(0, END)
            entry.insert(0, foldername)

    def run_script():
        """
        Start the analysis by calling the updated run_cats function.
        Note: Stop functionality has been removed because run_cats is not run as an external subprocess.
        """
        # Disable the Run button to avoid multiple clicks
        run_button.config(state='disabled')
        
        # Clear the output text box and reset progress bar
        output_text.configure(state='normal')
        output_text.delete('1.0', END)
        output_text.configure(state='disabled')
        progress_bar["value"] = 0

        # Update status
        status_label.config(text="Running...", foreground="red")
        status_label.grid()

        fasta_selection = fasta_combobox.get().strip()
        if fasta_selection == "CUSTOM":
            fasta = fasta_entry.get().strip()
        else:
            fasta = fasta_selection

        seq1 = seq1_entry.get()
        seq2 = seq2_entry.get()
        gtf = gtf_entry.get()

        # Build output path from folder, filename, and extension
        output_folder = output_folder_entry.get().strip()
        file_name = output_filename_entry.get().strip()
        extension = extension_var.get().strip()
        if not output_folder:
            messagebox.showerror("Error", "Please select an output folder.")
            status_label.config(text="", foreground="red")
            status_label.grid_remove()
            run_button.config(state='normal')
            return
        if not file_name:
            messagebox.showerror("Error", "Please specify a file name (without extension).")
            status_label.config(text="", foreground="red")
            status_label.grid_remove()
            run_button.config(state='normal')
            return
        if extension not in ["csv", "tsv"]:
            messagebox.showerror("Error", "Please select a valid extension: csv, tsv.")
            status_label.config(text="", foreground="red")
            status_label.grid_remove()
            run_button.config(state='normal')
            return

        final_output = os.path.join(output_folder, f"{file_name}.{extension}")

        # Get numeric parameters and convert to int if provided
        try:
            window_size_val = int(window_size_entry.get() or 5)
        except ValueError:
            messagebox.showerror("Error", "Window Size must be an integer.")
            run_button.config(state='normal')
            return

        try:
            num_bases_val = int(num_bases_entry.get() or 25)
        except ValueError:
            messagebox.showerror("Error", "Num Bases must be an integer.")
            run_button.config(state='normal')
            return

        variant_window_str = variant_window_entry.get().strip()
        if variant_window_str:
            try:
                variant_window_val = int(variant_window_str)
            except ValueError:
                messagebox.showerror("Error", "Variant Window must be an integer.")
                run_button.config(state='normal')
                return
        else:
            variant_window_val = None  # Will default to num_bases inside run_cats

        pathogenicity = pathogenicity_var.get()
        snv = snv_var.get()
        gene_list = gene_list_entry.get()

        def update_progress(current, total):
            progress_percent = int((current / total) * 100)
            progress_bar["value"] = progress_percent

        def progress_callback(current, total):
            root.after(0, update_progress, current, total)

        def run_in_thread():
            try:
                run_cats(
                    fasta_file=fasta,
                    seq1=seq1,
                    seq2=seq2,
                    gtf_file=gtf,
                    output=final_output,
                    window_size=window_size_val,
                    num_bases=num_bases_val,
                    pathogenicity=pathogenicity,
                    snv=snv,
                    gene_list=gene_list,
                    variant_window=variant_window_val,
                    progress_callback=progress_callback
                )
                messagebox.showinfo("Success", "Script ran successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Script failed: {e}")
            finally:
                status_label.config(text="", foreground="red")
                status_label.grid_remove()
                run_button.config(state='normal')
                root.after(0, update_progress, 0, 1)  # reset the progress bar

        threading.Thread(target=run_in_thread).start()

    def update_fields(event):
        """
        Automatically update or hide fields based on the FASTA dropdown selection.
        """
        selection = fasta_combobox.get().strip()
        if selection == "HUMAN":
            fasta_entry.grid_remove()
            browse_button.grid_remove()
            fasta_var.set("human")
            gtf_entry.delete(0, END)
            gtf_entry.insert(0, "../db/human/gencode.v47.annotation.gtf.gz")
        elif selection == "HUMAN Protein coding":
            fasta_entry.grid_remove()
            browse_button.grid_remove()
            fasta_var.set("human_pc")
            gtf_entry.delete(0, END)
            gtf_entry.insert(0, "../db/human/gencode.v47.annotation.gtf.gz")
        elif selection == "MOUSE":
            fasta_entry.grid_remove()
            browse_button.grid_remove()
            fasta_var.set("mouse")
            gtf_entry.delete(0, END)
            gtf_entry.insert(0, "../db/mouse/gencode.vM36.annotation.gtf.gz")
        elif selection == "MOUSE Protein coding":
            fasta_entry.grid_remove()
            browse_button.grid_remove()
            fasta_var.set("mouse_pc")
            gtf_entry.delete(0, END)
            gtf_entry.insert(0, "../db/mouse/gencode.vM36.annotation.gtf.gz")
        elif selection == "CUSTOM":
            fasta_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
            browse_button.grid(row=0, column=2, padx=5, pady=5)

    def update_pathogenicity(*args):
        """
        Ensure 'Pathogenic' is checked if 'SNVs' is checked.
        """
        if snv_var.get():
            pathogenicity_var.set(True)

    root = ttk.Window(themename="cosmo", iconphoto=ICON_FILE)
    root.title("Comparing Cas Activities by Target Superimposition (CATS)")

    if platform.system() == "Windows":
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)

    root.tk.call('tk', 'scaling', 1.4)

    root.geometry("1100x950")
    root.minsize(1100, 950)

    root.grid_rowconfigure(0, weight=0)
    root.grid_rowconfigure(1, weight=1)
    for c in range(4):
        root.grid_columnconfigure(c, weight=1)

    FONT = "Helvetica"
    font_large = (FONT, 17)
    font_medium = (FONT, 15)

    s = ttk.Style()
    s.configure('general.TButton', font=(FONT, 13))

    s2 = ttk.Style()
    s2.configure('run.TButton', font=(FONT, 13), background='green', foreground='white')
    s2.map(
        'run.TButton',
        background=[('active', 'lightgreen'), ('!active', 'green')],
        foreground=[('active', 'white'), ('!active', 'white')]
    )

    s4 = ttk.Style()
    s4.configure('general.success.TCheckbutton', font=font_large, indicatorsize=20)

    s5 = ttk.Style()
    s5.configure('general.TCombobox', font=(FONT, 17, "bold"), background='white')

    s.configure('Banner.TLabel', font=(FONT, 24, "bold"), foreground='blue')

    root.option_add('*TCombobox*Listbox.font', font_medium)

    banner_frame = ttk.Frame(root, padding=(10, 10, 10, 10))
    banner_frame.grid(row=0, column=0, columnspan=4, sticky="nswe")
    banner_frame.config(style='Banner.TFrame')

    original_image = Image.open(ICON_FILE)
    scale_factor = 0.10
    width, height = original_image.size
    resized_image = original_image.resize((int(width * scale_factor), int(height * scale_factor)),
                                           Image.Resampling.LANCZOS)
    cat_image = ImageTk.PhotoImage(resized_image)

    banner_label = ttk.Label(
        banner_frame,
        text="    Welcome to CATS",
        style='Banner.TLabel',
        image=cat_image,
        compound='left',
        padding=5
    )
    banner_label.image = cat_image
    banner_label.pack(fill="x")

    notebook = ttk.Notebook(root, bootstyle="secondary")
    notebook.grid(row=1, column=0, columnspan=4, sticky="nsew", padx=10, pady=10)

    input_frame = ttk.Frame(notebook)
    logging_frame = ttk.Frame(notebook)
    docs_frame = ttk.Frame(notebook)

    notebook.add(input_frame, text=" CATS Input ")
    notebook.add(logging_frame, text=" Logging ")
    notebook.add(docs_frame, text=" Docs ")

    docs_title = ttk.Label(
        docs_frame,
        text="Comparing Cas Activities by Target Superimposition (CATS)",
        font=(FONT, 20, "bold")
    )
    docs_title.grid(row=0, column=0, padx=10, pady=10, sticky="nw")

    docs_text = ttk.Text(docs_frame, wrap="word", font=font_medium)
    docs_text.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
    docs_text.tag_configure("title", font=("Helvetica", 16, "bold", "underline"), spacing3=10)
    docs_text.tag_configure("subtitle", font=("Helvetica", 14, "bold"), foreground="#333")
    docs_text.tag_configure("body", font=("Helvetica", 12), spacing3=6)
    docs_text.tag_configure("bulleted", font=("Helvetica", 12), lmargin1=40, lmargin2=40)

    docs_text.insert("end", "What is CATS?\n", "title")
    docs_text.insert("end",
        "CATS is a bioinformatic tool designed to identify genomic regions near "
        "two distinct sequences of interest.\n\n", "body"
    )

    docs_text.insert("end", "CATS Input tab Documentation:\n", "title")
    docs_text.insert("end", 
        "The CATS Input Tab is where you configure the parameters "
        "and settings for running CATS. Below are the main parameters.\n", "body"
    )

    docs_text.insert("end", "Genome Inputs:\n", "subtitle")
    docs_text.insert("end", 
        "-  FASTA Dropdown: Choose a genome reference from the"
        " following options: 'HUMAN', 'HUMAN Protein coding', 'MOUSE', "
        "'MOUSE Protein coding', 'CUSTOM'.\n", "bulleted"
    )
    docs_text.insert("end", 
        "The HUMAN and MOUSE choices will automatically retrieve"
        " the corresponding GENCODE genomes and annotations. When selecting 'CUSTOM',"
        " a text box and 'Browse' button appear, allowing the user to provide a"
        " custom FASTA file path.\n", "bulleted"
    )
    docs_text.insert("end", 
        "-  GTF Entry: Automatically filled based on the chosen genome. "
        "You can also manually browse for a GTF file. It is not necessary"
        " for a custom FASTA file.\n\n", "bulleted"
    )

    docs_text.insert("end", "Sequences:\n", "subtitle")
    docs_text.insert("end",
        "-  Sequence 1: Enter the first (or unique) PAM sequence"
        " for the analysis. Required field.\n", "bulleted"
    )
    docs_text.insert("end", 
        "-  Sequence 2 (optional): Enter a secondary sequence for"
        " automatic detection of overlapping PAM sequences.\n\n", "bulleted"
    )

    docs_text.insert("end", "Output:\n", "subtitle")
    docs_text.insert("end",
        "-  Output Folder: Specify the directory where results will be saved. "
        "Use the 'Browse' button to select a folder. Required field.\n", "bulleted"
    )
    docs_text.insert("end",
        "-  Filename: Enter the base name of the output file (without extension). "
        "Required field.\n", "bulleted"
    )
    docs_text.insert("end", 
        "-  Extension: Choose between CSV or TSV formats for "
        "the output file.\n\n", "bulleted"
    )

    docs_text.insert("end", "Parsing parameters:\n", "subtitle")
    docs_text.insert("end",
        "-  Window Size (default 5): Set the size of the window "
        "around the sequences (for double-sequence mode).\n", "bulleted"
    )
    docs_text.insert("end",
        "-  Num Bases (default 25): Specify the number of bases"
        " preceding and succeeding each sequence to include in the output.\n", "bulleted"
    )
    docs_text.insert("end",
        "-  Gene List: Optionally specify a file containing a list of gene names "
        "to include in the analysis. You can also enter gene names directly, separated"
        " by semicolons (e.g., HBB;HTT).\n\n", "bulleted"
    )

    docs_text.insert("end", "Pathogenic variants:\n", "subtitle")
    docs_text.insert("end",
        "-  Pathogenic: Check to include only sequences containing potentially"
        " pathogenic variants, as identified by ClinVar.\n", "bulleted"
    )
    docs_text.insert("end",
        "-  SNVs: Check to include only sequences associated with"
        " single nucleotide variants (SNVs) from ClinVar. Enabling this also"
        " checks 'Pathogenic'.\n", "bulleted"
    )
    docs_text.insert("end",
        "-  Variant Window: Specify the maximum distance between a"
        " mutation and the found PAM sequence. If left blank, it will "
        "use the number of bases specified earlier. A distance of 0 corresponds"
        " to a mutation inside the PAM sequence.\n\n", "bulleted"
    )

    docs_text.insert("end",
        "After configuring the inputs, click the Run button. "
        "Monitor the process in the Logging tab.\n\n", "body"
    )
    docs_text.insert("end", "GiHub repository:  ", "body")
    docs_text.insert("end", "https://github.com/Physics4MedicineLab/CATS\n", "body")

    docs_text.configure(state="disabled")

    for i in range(12):
        input_frame.grid_rowconfigure(i, weight=0)
        input_frame.grid_columnconfigure(i, weight=1)
    logging_frame.grid_rowconfigure(0, weight=1)
    logging_frame.grid_columnconfigure(0, weight=1)

    genome_frame = ttk.Labelframe(input_frame, text="  Genome Inputs  ", bootstyle="primary")
    genome_frame.grid(row=0, column=0, columnspan=4, sticky="nsew", padx=5, pady=5)

    genome_frame.grid_columnconfigure(0, weight=0, minsize=180)
    genome_frame.grid_columnconfigure(1, weight=1)
    genome_frame.grid_columnconfigure(2, weight=0, minsize=80)

    fasta_var = ttk.StringVar()
    gtf_var = ttk.StringVar()

    fasta_combobox = ttk.Combobox(
        genome_frame,
        textvariable=fasta_var,
        font=font_large,
        style='general.TCombobox',
        values=["HUMAN", "HUMAN Protein coding", "MOUSE", "MOUSE Protein coding", "CUSTOM"],
        state="readonly"
    )
    fasta_combobox.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
    fasta_combobox.bind("<<ComboboxSelected>>", update_fields)

    fasta_entry = ttk.Entry(genome_frame, font=font_large, width=50)
    browse_button = ttk.Button(
        genome_frame,
        text="Browse",
        command=lambda: browse_file(fasta_entry),
        style='general.TButton'
    )

    gtf_label = ttk.Label(genome_frame, text="GTF:", font=font_large)
    gtf_label.grid(row=1, column=0, sticky="w", padx=5, pady=5)
    gtf_entry = ttk.Entry(genome_frame, textvariable=gtf_var, font=font_large, width=50)
    gtf_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=5)
    gtf_browse_button = ttk.Button(genome_frame, text="Browse", command=lambda: browse_file(gtf_entry), style='general.TButton')
    gtf_browse_button.grid(row=1, column=2, padx=5, pady=5)

    seq_frame = ttk.Labelframe(input_frame, text="  Sequences  ", bootstyle="primary")
    seq_frame.grid(row=1, column=0, columnspan=4, sticky="nsew", padx=5, pady=5)

    seq_frame.grid_columnconfigure(0, weight=0, minsize=180)
    seq_frame.grid_columnconfigure(1, weight=1)

    seq1_var = ttk.StringVar()
    seq2_var = ttk.StringVar()

    def create_labeled_entry(master, text, row, column, var):
        label = ttk.Label(master, text=text, font=font_large)
        label.grid(row=row, column=column, sticky="w", padx=5, pady=5)
        entry = ttk.Entry(master, textvariable=var, font=font_large, width=50)
        entry.grid(row=row, column=column+1, sticky="ew", padx=5, pady=5)
        return entry

    seq1_entry = create_labeled_entry(seq_frame, "Sequence 1:", 0, 0, seq1_var)
    seq2_entry = create_labeled_entry(seq_frame, "Sequence 2 [optional]:   ", 1, 0, seq2_var)

    output_frame = ttk.Labelframe(input_frame, text="  Output  ", bootstyle="primary")
    output_frame.grid(row=2, column=0, columnspan=4, sticky="nsew", padx=5, pady=5)

    output_frame.grid_columnconfigure(0, weight=0, minsize=180)
    output_frame.grid_columnconfigure(1, weight=1)
    output_frame.grid_columnconfigure(2, weight=0, minsize=80)
    output_frame.grid_columnconfigure(3, weight=0, minsize=80)

    output_folder_var = ttk.StringVar()
    output_folder_label = ttk.Label(output_frame, text="Output Folder:               ", font=font_large)
    output_folder_label.grid(row=0, column=0, sticky="w", padx=5, pady=5)

    output_folder_entry = ttk.Entry(output_frame, textvariable=output_folder_var, font=font_large, width=40)
    output_folder_entry.grid(row=0, column=1, sticky="ew", padx=5, pady=5, columnspan=2)

    output_folder_browse = ttk.Button(output_frame, text="Browse", command=lambda: browse_folder(output_folder_entry), style='general.TButton')
    output_folder_browse.grid(row=0, column=3, padx=5, pady=5)

    output_filename_var = ttk.StringVar()
    output_filename_label = ttk.Label(output_frame, text="Filename:", font=font_large)
    output_filename_label.grid(row=1, column=0, sticky="w", padx=5, pady=5)

    output_filename_entry = ttk.Entry(output_frame, textvariable=output_filename_var, font=font_large, width=25)
    output_filename_entry.grid(row=1, column=1, sticky="ew", padx=5, pady=5)

    extension_var = ttk.StringVar()
    extension_label = ttk.Label(output_frame, text="Extension:", font=font_large)
    extension_label.grid(row=1, column=2, sticky="w", padx=5, pady=5)
    extension_combobox = ttk.Combobox(output_frame, textvariable=extension_var, values=["csv", "tsv"], font=font_large, width=8, state="readonly", style='general.TCombobox')
    extension_combobox.set("csv")
    extension_combobox.grid(row=1, column=3, sticky="w", padx=5, pady=5)

    param_frame = ttk.Labelframe(input_frame, text="  Parsing parameters  ", bootstyle="primary")
    param_frame.grid(row=3, column=0, columnspan=4, sticky="nsew", padx=5, pady=5)
    param_frame.grid_columnconfigure(0, weight=0, minsize=180)
    param_frame.grid_columnconfigure(1, weight=1)
    param_frame.grid_columnconfigure(2, weight=0, minsize=80)

    window_size_var = ttk.StringVar()
    num_bases_var = ttk.StringVar()

    window_size_entry = create_labeled_entry(param_frame, "Window Size (default 5):", 2, 0, window_size_var)
    num_bases_entry = create_labeled_entry(param_frame, "Num Bases (default 25):", 3, 0, num_bases_var)

    gene_list_var = ttk.StringVar()
    gene_list_label = ttk.Label(param_frame, text="Gene List:", font=font_large)
    gene_list_label.grid(row=4, column=0, sticky="w", padx=5, pady=5)
    gene_list_entry = ttk.Entry(param_frame, textvariable=gene_list_var, font=font_large, width=50)
    gene_list_entry.grid(row=4, column=1, sticky="ew", padx=5, pady=5)
    gene_list_browse = ttk.Button(param_frame, text="Browse", command=lambda: browse_file(gene_list_entry), style='general.TButton')
    gene_list_browse.grid(row=4, column=2, padx=5, pady=5)

    flag_frame = ttk.Labelframe(input_frame, text="  Pathogenic variants - Human only  ", bootstyle="primary")
    flag_frame.grid(row=4, column=0, columnspan=4, sticky="nsew", padx=5, pady=5)
    flag_frame.grid_columnconfigure(0, weight=0, minsize=180)
    flag_frame.grid_columnconfigure(1, weight=1)

    pathogenicity_var = ttk.BooleanVar()
    snv_var = ttk.BooleanVar()

    # Ensure 'Pathogenic' is checked if 'SNVs' is checked
    snv_var.trace_add("write", update_pathogenicity)
    pathogenicity_checkbox = ttk.Checkbutton(flag_frame, text="Pathogenic", variable=pathogenicity_var, style="general.success.TCheckbutton")
    pathogenicity_checkbox.grid(row=0, column=0, padx=20, pady=5, sticky="w")
    snv_checkbox = ttk.Checkbutton(flag_frame, text="SNVs", variable=snv_var, style="general.success.TCheckbutton")
    snv_checkbox.grid(row=0, column=1, padx=20, pady=5, sticky="w")

    variant_window_var = ttk.StringVar()
    variant_window_entry = create_labeled_entry(flag_frame, "Variant Window:            ", 1, 0, variant_window_var)

    button_frame = ttk.Frame(input_frame)
    button_frame.grid(row=6, column=0, columnspan=4, sticky="nsew", padx=5, pady=20)
    run_button = ttk.Button(button_frame, text="Run", command=run_script, style='run.TButton')
    run_button.grid(row=0, column=0, padx=10)

    status_label = ttk.Label(button_frame, text="", font=font_large, foreground="red")
    status_label.grid(row=0, column=1, padx=10, sticky="nsew")
    status_label.grid_remove()

    output_text = ttk.Text(logging_frame, wrap='word', state='disabled', font=font_medium)
    output_text.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

    scrollbar = ttk.Scrollbar(logging_frame, orient='vertical', command=output_text.yview)
    scrollbar.grid(row=0, column=1, sticky='nsew')
    output_text.config(yscrollcommand=scrollbar.set)

    progress_bar = ttk.Progressbar(logging_frame, orient="horizontal", mode="determinate")
    progress_bar.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")
    progress_bar["maximum"] = 100

    output_text.tag_configure("info", foreground="orange")
    output_text.tag_configure("warning", foreground="red")
    output_text.tag_configure("timestamp", foreground="green")
    output_text.tag_configure("message", foreground="black")
    output_text.configure(foreground="black")

    sys.stdout = TextRedirector(output_text)
    sys.stderr = TextRedirector(output_text)

    root.mainloop()

if __name__ == "__main__":
    main()
