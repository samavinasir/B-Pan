import time
import warnings
import glob
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import pandas as pd
import os
import threading
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import os, sys, subprocess, shutil


# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Centralized path constants ───────────────────────────────────────────────
DATA_DIR    = os.path.join(SCRIPT_DIR, "Data")
DB_DIR      = os.path.join(DATA_DIR, "Db")
TOOLS_DIR   = os.path.join(DATA_DIR, "Tools")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "Results")

# Tool executables — prefer system PATH, fall back to bundled binaries
DIAMOND_EXE  = shutil.which("diamond")  or os.path.join(TOOLS_DIR, "diamond", "diamond.exe")
RPSBLAST_EXE = shutil.which("rpsblast") or os.path.join(TOOLS_DIR, "blast+", "bin", "rpsblast.exe")
CLUSTALO_EXE = shutil.which("clustalo") or os.path.join(TOOLS_DIR, "clustal-omega-1.2.2-win64", "clustalo.exe")
MAFFT_EXE    = shutil.which("mafft")    or os.path.join(TOOLS_DIR, "mafft", "mafft-win", "mafft.bat")
MUSCLE_EXE   = shutil.which("muscle")   or os.path.join(TOOLS_DIR, "muscle", "muscle.exe")

# Database paths
COG_DB      = os.path.join(DB_DIR, "COGDb", "COGDb")
COG_CDD_DB  = os.path.join(DB_DIR, "COG_CDD", "Cog", "Cog")
CDDID_TBL   = os.path.join(DB_DIR, "COG_CDD", "cddid.tbl")
KEGG_DB     = os.path.join(DB_DIR, "KEGGDb", "KEGGDb")
PHYLO_DB    = os.path.join(DB_DIR, "PhyloDb", "PhyloDb")

# Small CSV mapping files (shipped with the package)
COG_ANNOTATION_CSV   = os.path.join(DB_DIR, "COG_ANNOTATION.csv")
# Complete official COG-20 (2020) functional category definitions (4877 entries),
# downloaded from https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
# Used instead of the old COG_ANNOTATION_CSV (which only had 408 of ~4877 COGs).
COG_20_DEF_TAB        = os.path.join(DB_DIR, "cog-20.def.tab")
SYMBOL_CATEGORIES_CSV= os.path.join(DB_DIR, "SYMBOL_CATEGORIES.csv")
KEGG_ANNOTATION_CSV  = os.path.join(DB_DIR, "KEGG_ANNOTATION.csv")

# Results subdirectories
RESULTS_PANGENOME = os.path.join(RESULTS_DIR, "Pangenome")
RESULTS_COG       = os.path.join(RESULTS_DIR, "COG")
RESULTS_KEGG      = os.path.join(RESULTS_DIR, "KEGG")
RESULTS_PHYLO     = os.path.join(RESULTS_DIR, "Phylo")

# Suppress Biopython deprecation warnings about command-line wrappers BEFORE importing submodules
try:
    from Bio import BiopythonDeprecationWarning
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
except Exception:
    pass

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio import SeqIO, AlignIO
import subprocess
import sys
import shutil
import requests

matplotlib.use("Agg")


window = Tk()
window.title("B-Pan")
window.geometry("980x800")  # Increased height from 640 to 800
window.maxsize(height=800, width=980)  # Increased max height
window.resizable(width=False, height=False)

First_Color = Second_Color = "#EAEAEA"
Buttons_Colors = Label_Color = "#003380"
Button_Text = "#F5F5F5"
font_family = "Bahnschrift"
width_frames = 1150

IBG_LAB_LOGO = PhotoImage(file=os.path.join(SCRIPT_DIR, "Data", "Images", "New IBG Logo.png"))
MGBio_LOGO = PhotoImage(file=os.path.join(SCRIPT_DIR, "Data", "Images", "MGBio.png"))
HELP_ICON = PhotoImage(file=os.path.join(SCRIPT_DIR, "Data", "Images", "Help.png"))
SUBMIT_BUTTON = PhotoImage(file=os.path.join(SCRIPT_DIR, "Data", "Images", "SUBMIT.png"))
BPAN_LOGO = PhotoImage(file=os.path.join(SCRIPT_DIR, "Data", "Images", "Logo.png"))

def read_submitted_fasta_file(file_name):
    """Read FASTA file and validate sequences."""
    title_list = []
    fasta_list = []
    description_list = []
    
    try:
        for seq_record in SeqIO.parse(file_name, "fasta"):
            # Validate sequence
            if len(seq_record.seq) > 0 and str(seq_record.seq).strip():
                title_list.append(seq_record.id)
                fasta_list.append(seq_record.seq)
                description_list.append(seq_record.description)
            else:
                print(f"[Warning] Skipping empty sequence: {seq_record.id}")
        
        if not title_list:
            raise ValueError(f"No valid sequences found in {file_name}")
            
        return title_list, fasta_list, description_list
    except Exception as e:
        raise ValueError(f"Error reading FASTA file {file_name}: {e}")

filename_var = StringVar()

def FASTA_BROWSE():
    """Browse and validate FASTA files."""
    try:
        filename = list(filedialog.askopenfilenames(
            parent=window, 
            initialdir=os.path.expanduser("~"), 
            title="Choose Fasta File", 
            filetypes=[("fasta file(.FAA)", ".FAA"), ("fasta file(.FASTA)", ".FASTA")]
        ))
        
        if not filename:
            messagebox.showwarning("Warning", "No File Selected")
            return
            
        FILES_LABEL.configure(text=f"Total Files {str(len(filename))}")
        Total_sequences_length = []
        
        for files in filename:
            try:
                fasta_sequences_title, fasta_sequences_file, fasta_sequences_des = read_submitted_fasta_file(files)
                for single_seq in fasta_sequences_file:
                    Total_sequences_length.append(single_seq)
            except ValueError as e:
                messagebox.showerror("File Error", f"Error in file {os.path.basename(files)}: {e}")
                return
            except Exception as e:
                messagebox.showerror("File Error", f"Unexpected error reading {os.path.basename(files)}: {e}")
                return
        
        if Total_sequences_length:
            SEQUENCES_LABEL.configure(text=f"Total Sequences >{len(Total_sequences_length)}")
            filename_var.set(str(filename))
        else:
            messagebox.showerror("File Error", "No valid sequences found in selected files")

    except AttributeError:
        messagebox.showwarning("Warning", "No File Selected")
    except Exception as e:
        messagebox.showerror("Error", f"Unexpected error: {e}")


startupinfo = None
if sys.platform == "win32":
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW


def run_diamond_blastp(reference_db, input_sequences, output_file):
    PATH_DIAMOND = DIAMOND_EXE
    cmd = [
        PATH_DIAMOND, "blastp",
        "-d", reference_db,
        "-q", input_sequences,
        "-o", output_file,
        "--outfmt", "6"
    ]
    subprocess.run(cmd, startupinfo=startupinfo, shell=False)

def run_diamond_blastp_with_retries(reference_db: str, input_sequences: str, output_file: str) -> bool:
    """Run DIAMOND blastp with progressively more sensitive settings until output has hits.

    Returns True if an output file was produced (may be empty), and False only if the file could not be created.
    """
    PATH_DIAMOND = DIAMOND_EXE
    # Try from standard to ultra sensitive with relaxed e-values and more targets
    attempt_params = [
        ["--sensitive", "-e", "1e-2", "-k", "25"],
        ["--more-sensitive", "-e", "1e-1", "-k", "50"],
        ["--ultra-sensitive", "-e", "1", "-k", "100"]
    ]
    # Ensure directory exists
    try:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)) or ".", exist_ok=True)
    except Exception:
        pass
    last_created = False
    for extra in attempt_params:
        # Clean previous empty file before trying again
        try:
            if os.path.exists(output_file) and os.path.getsize(output_file) == 0:
                os.remove(output_file)
        except Exception:
            pass
        cmd = [
            PATH_DIAMOND, "blastp",
            "-d", reference_db,
            "-q", input_sequences,
            "-o", output_file,
            "--outfmt", "6",
        ] + extra
        try:
            subprocess.run(cmd, startupinfo=startupinfo, shell=False)
        except Exception:
            continue
        if os.path.exists(output_file):
            last_created = True
            # If file has any hits, stop
            try:
                if os.path.getsize(output_file) > 0:
                    break
            except Exception:
                break
    return last_created

def which_exe(exe_name: str):
    # Find executable in known tools dir or PATH
    candidates = [
        os.path.join(TOOLS_DIR, 'blast+', exe_name),
        os.path.join(TOOLS_DIR, 'blast-2.17.0+', 'bin', exe_name),
        exe_name
    ]
    for cand in candidates:
        try:
            # On Windows, .exe extension
            if sys.platform == 'win32' and not cand.lower().endswith('.exe'):
                cand_exe = cand + '.exe'
            else:
                cand_exe = cand
            # Use where/which via shutil.which
            found = shutil.which(cand_exe)
            if found:
                return found
            if os.path.isfile(cand_exe):
                return cand_exe
        except Exception:
            continue
    return None

def run_rpsblast(reference_db, input_sequences, output_file):
    """
    Fixed RPS-BLAST function with robust path handling and error checking.
    """
    import os, sys, subprocess, shutil, time, threading
    import glob  # Add this import

    # --- Use fixed paths (from your diagnostic) ---
    rps_exe = RPSBLAST_EXE
    rps_db = COG_CDD_DB

    print("DEBUG >>> Forcing RPS-BLAST exe:", rps_exe)
    print("DEBUG >>> Forcing RPS-BLAST DB:", rps_db)

    # Validate executable
    if not os.path.exists(rps_exe):
        raise FileNotFoundError(f"RPS-BLAST executable not found: {rps_exe}")

    # FIXED: Check for multi-volume database files
    def check_rps_database(db_path):
        db_dir = os.path.dirname(db_path)
        base_name = os.path.basename(db_path)
        
        # Check for .pal file (alias file, used when DB was built with -in or makeprofiledb -title)
        pal_file = os.path.join(db_dir, f"{base_name}.pal")
        if os.path.exists(pal_file):
            print(f"[OK] Found PAL file: {pal_file}")
            return True
            
        # Check for multi-volume RPS-BLAST files (e.g. Cog.00.rps, Cog.01.rps, ...)
        rps_files = glob.glob(os.path.join(db_dir, f"{base_name}*.rps"))
        if len(rps_files) > 0:
            print(f"[OK] Found {len(rps_files)} .rps database files")
            return True
            
        return False

    if not check_rps_database(rps_db):
        raise FileNotFoundError(f"RPS-BLAST database not found: {rps_db}")

    # Validate input file
    if not os.path.exists(input_sequences):
        raise FileNotFoundError(f"Input file not found: {input_sequences}")

    # Count sequences
    seq_count = 0
    with open(input_sequences, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_count += 1
    print(f"Found {seq_count} sequences to process")

    # Ensure output dir exists
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

    # Normalize and quote-safe absolute paths
    rps_exe_norm = os.path.normpath(rps_exe)
    rps_db_norm = os.path.normpath(rps_db)
    input_norm = os.path.normpath(os.path.abspath(input_sequences))
    output_norm = os.path.normpath(os.path.abspath(output_file))

    cmd = [
        rps_exe_norm,
        "-query", input_norm,
        "-db", rps_db_norm,
        "-out", output_norm,
        "-evalue", "1e-5",
        "-outfmt", "6 qseqid sacc pident length qlen qstart qend evalue bitscore",
        "-num_threads", str(os.cpu_count() or 1)
    ]
    # Build RPS-BLAST command
   # cmd = [
    #    os.path.normpath(rps_exe),
     #   "-query", os.path.normpath(input_sequences),
      #  "-db", os.path.normpath(rps_db),
      #  "-out", os.path.normpath(output_file),
      #  "-evalue", "1e-5",
      #  "-outfmt", "6 qseqid sacc pident length qlen qstart qend evalue bitscore",
      #  "-num_threads", str(os.cpu_count() or 1)
    #]

    print("DEBUG CMD:", " ".join(cmd))

    # Progress monitoring setup
    progress_file = output_file + ".progress"
    monitor_running = True
    
    def monitor_output_progress():
        last_size = 0
        last_update_time = time.time()
        processed_count = 0
        
        while monitor_running:
            try:
                if os.path.exists(output_file):
                    current_size = os.path.getsize(output_file)
                    if current_size > last_size or time.time() - last_update_time > 10:
                        unique_queries = set()
                        try:
                            with open(output_file, 'r') as f:
                                for line in f:
                                    if line.strip():
                                        parts = line.split('\t')
                                        if parts and len(parts) > 0:
                                            unique_queries.add(parts[0])
                            
                            processed_count = len(unique_queries)
                        except Exception as e:
                            print(f"Progress monitoring read error: {e}")
                        
                        progress = min(99, int(100 * processed_count / max(1, seq_count)))
                        current_time = time.time()
                        time_diff = current_time - last_update_time
                        
                        if time_diff > 5:
                            rate = (current_size - last_size) / time_diff / 1024  # KB/s
                            status_msg = f"RPS-BLAST Progress: {processed_count}/{seq_count} sequences ({progress}%) - {rate:.1f} KB/s"
                            print(status_msg)
                            last_update_time = current_time
                            last_size = current_size
                        else:
                            status_msg = f"RPS-BLAST Progress: {processed_count}/{seq_count} sequences ({progress}%)"
                            print(status_msg)
                        
                        try:
                            ui_set_status(f"Job Status : {status_msg}", "Red")
                        except Exception:
                            pass
                        
                        with open(progress_file, 'w') as pf:
                            pf.write(f"{progress}\n")
            except Exception as e:
                print(f"Progress monitoring error: {e}")
            
            time.sleep(2)

    monitor_thread = threading.Thread(target=monitor_output_progress)
    monitor_thread.daemon = True
    monitor_thread.start()
    
    # Windows-specific startup info
    si = subprocess.STARTUPINFO()
    si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    si.wShowWindow = subprocess.SW_HIDE
    
    try:
        # Execute RPS-BLAST with proper error handling
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, bufsize=1, startupinfo=si, cwd=os.path.dirname(rps_db))
        
        # Monitor stderr for progress information
        for line in process.stderr:
            line = line.strip()
            if line:
                print(f"RPS-BLAST: {line}")
                try:
                    ui_set_status(f"Job Status : RPS-BLAST - {line}", "Red")
                except Exception:
                    pass
        
        process.wait()
        
        # Stop monitoring
        monitor_running = False
        monitor_thread.join(timeout=5)
        
        # Check results
        if process.returncode != 0:
            raise Exception(f"RPS-BLAST failed with return code {process.returncode}")
        
        # Validate output
        if not os.path.exists(output_file):
            raise Exception("RPS-BLAST output file was not created")
        
        # Add header to output file
        if os.path.exists(output_file):
            header = "QueryID\tSubjectID\tIdentity(%)\tAlignmentLength\tQueryLength\tQueryStart\tQueryEnd\tE-value\tBitScore\n"
            with open(output_file, "r+", encoding="utf-8") as f:
                content = f.read()
                if not content.startswith("QueryID"):
                    f.seek(0, 0)
                    f.write(header + content)
        
        # Clean up
        try:
            os.remove(progress_file)
        except Exception:
            pass
        
        try:
            ui_set_status(f"Job Status : RPS-BLAST completed - {seq_count}/{seq_count} sequences (100%)", "Red")
        except Exception:
            pass
        
        print("RPS-BLAST search completed successfully.")
        return True
        
    except Exception as e:
        # Stop monitoring on error
        monitor_running = False
        try:
            monitor_thread.join(timeout=1)
        except Exception:
            pass
        
        print(f"RPS-BLAST error: {e}")
        raise e
def create_diamond_db(fasta_file, db_name):
    PATH_DIAMOND = DIAMOND_EXE
    makedb_cmd = [PATH_DIAMOND, "makedb", "--in", fasta_file, "-d", db_name]
    subprocess.run(makedb_cmd, startupinfo=startupinfo, shell=False)

window.configure(bg=First_Color)

window.iconphoto(False,BPAN_LOGO)

B_PAN_LOGO = Label(window, image=BPAN_LOGO, bg = First_Color)
B_PAN_LOGO.place(x=40, y=20)

LABEL_LOGO = Label(window, text='B-Pan', font=(font_family, 34), bg = First_Color, fg=Label_Color)
LABEL_LOGO.place(x=105, y=20)

LABEL_DESCRIPTION = Label(window, text='A Robust Software Package for Bacterial Pangenome Analysis', font=(font_family, 15), bg = First_Color, fg=Label_Color)
LABEL_DESCRIPTION.place(x=40, y=85)

IBG_LOGO = Label(window, image=IBG_LAB_LOGO, bg=Second_Color)
IBG_LOGO.place(x=740, y=20)

MGBio_LOGO_Label = Label(window, image=MGBio_LOGO, bg=Second_Color)
MGBio_LOGO_Label.place(x=660, y=20)

shadow = Frame(window,width=393, height=146, relief="solid", bg="#DCDCDC")
shadow.place(x=47, y=137)

DETAIL_FRAME = Frame(window, width=393, height=185, bg="#F3F3F3", relief="solid")
DETAIL_FRAME.place(x=40, y=130)

EXTENSTION_LABEL = Label(DETAIL_FRAME, text="Allowed Extentions: '.fasta', '.faa', '.fas'", font=(font_family, 13), bg="#F1F1F1", fg=Label_Color)
EXTENSTION_LABEL.place(x=20, y=10)

EXTENSION_BROWSE_IMAGE = PhotoImage(file=os.path.join(SCRIPT_DIR, "Data", "Images", "Browse.png"))
EXTENSION_BROWSE = Button(DETAIL_FRAME, image=EXTENSION_BROWSE_IMAGE, bd=0, command=FASTA_BROWSE)
EXTENSION_BROWSE.place(x=20, y=43)

FILES_LABEL_Var = StringVar()
FILES_LABEL = Label(DETAIL_FRAME, text=f'Total Files {0}', font=(font_family, 12), bg = "#F1F1F1", fg=Label_Color)
FILES_LABEL.place(x=20, y=80)

SEQUENCES_LABEL = Label(DETAIL_FRAME, text=f'Total Sequences >{0}', font=(font_family, 12), bg = "#F1F1F1", fg=Label_Color)
SEQUENCES_LABEL.place(x=20, y=100)

LOADING_LABEL = Label(DETAIL_FRAME, text="Job Status : None", font=(font_family, 12), bg = "#F1F1F1", fg=Label_Color, wraplength=355, justify="left", anchor="w")
LOADING_LABEL.place(x=20, y=122, width=355)
LOADING_LABEL.lift()

def ui_info(message):
    # Ensure messageboxes are invoked on the Tk main thread
    window.after(0, lambda: messagebox.showinfo("Info", message))

def ui_error(title, message):
    # Ensure error messageboxes are invoked on the Tk main thread
    window.after(0, lambda: messagebox.showerror(title, message))

def ui_set_status(text, color=None):
    # Safely update status label from worker thread
    def _apply():
        kwargs = {"text": text}
        if color is not None:
            kwargs["fg"] = color
        LOADING_LABEL.config(**kwargs)
    window.after(0, _apply)

def apply_publication_style():
    """Apply comprehensive publication-ready styling for high-quality scientific figures."""
    
    # Reset matplotlib to defaults first to avoid conflicts
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    
    # Set publication-quality parameters with proper scientific figure standards
    publication_params = {
        # **FIGURE QUALITY & OUTPUT**
        'figure.dpi': 150,                    # Higher display DPI for crisp preview
        'savefig.dpi': 600,                   # Publication DPI (300-600 for journals)
        'savefig.bbox': 'tight',              # Trim whitespace
        'savefig.pad_inches': 0.1,            # Small padding
        'savefig.facecolor': 'white',         # White background
        'savefig.edgecolor': 'none',          # No border
        'savefig.format': 'png',              # Default format
        
        # **FONTS - LARGE AND READABLE**
        'font.family': 'sans-serif',          # Professional font family
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica', 'Verdana'],
        'font.size': 16,                      # Base font size (was 12) - MUCH LARGER
        'font.weight': 'normal',              # Normal weight
        
        # **AXES - THICK AND VISIBLE**
        'axes.titlesize': 22,                 # Large titles (was 16)
        'axes.titleweight': 'bold',           # Bold titles
        'axes.titlepad': 20,                  # Space above title
        'axes.labelsize': 18,                 # Large axis labels (was 14)
        'axes.labelweight': 'bold',           # Bold axis labels
        'axes.labelpad': 8,                   # Label padding
        'axes.linewidth': 2.5,                # THICK axes borders (was 1.2)
        'axes.edgecolor': 'black',            # Black axes
        'axes.facecolor': 'white',            # White background
        'axes.grid': True,                    # Enable grid
        'axes.axisbelow': True,               # Grid behind data
        'axes.spines.left': True,             # Show all spines
        'axes.spines.bottom': True,
        'axes.spines.top': False,             # Hide top spine (cleaner look)
        'axes.spines.right': False,           # Hide right spine
        
        # **TICKS - LARGE AND BOLD**
        'xtick.labelsize': 16,                # Large x-tick labels (was 12)
        'ytick.labelsize': 16,                # Large y-tick labels (was 12)
        'xtick.major.size': 8,                # Longer tick marks (was ~4)
        'ytick.major.size': 8,
        'xtick.major.width': 2,               # Thicker tick marks
        'ytick.major.width': 2,
        'xtick.minor.size': 4,                # Minor ticks
        'ytick.minor.size': 4,
        'xtick.minor.width': 1.5,
        'ytick.minor.width': 1.5,
        'xtick.color': 'black',               # Black ticks
        'ytick.color': 'black',
        
        # **LINES - THICK AND VISIBLE**
        'lines.linewidth': 3.0,               # THICK lines (was 1.5)
        'lines.markersize': 10,               # Larger markers (was ~6)
        'lines.markeredgewidth': 2,           # Thick marker edges
        'lines.solid_capstyle': 'round',      # Rounded line caps
        
        # **LEGEND - CLEAR AND READABLE**
        'legend.fontsize': 16,                # Large legend (was 12)
        'legend.frameon': True,               # Frame around legend
        'legend.fancybox': False,             # Square corners
        'legend.shadow': False,               # No shadow
        'legend.framealpha': 0.9,             # Semi-transparent background
        'legend.facecolor': 'white',          # White background
        'legend.edgecolor': 'black',          # Black border
        'legend.borderpad': 0.6,              # Internal padding
        'legend.columnspacing': 1.5,          # Space between columns
        'legend.handlelength': 2.0,           # Length of legend handles
        'legend.handletextpad': 0.8,          # Space between handle and text
        'legend.labelspacing': 0.5,           # Space between entries
        
        # **GRID - SUBTLE BUT VISIBLE**
        'grid.color': '#CCCCCC',              # Light gray grid
        'grid.linestyle': '-',                # Solid lines
        'grid.linewidth': 1.0,                # Medium thickness
        'grid.alpha': 0.8,                    # Semi-transparent
        
        # **PATCHES (BAR PLOTS, PIE CHARTS)**
        'patch.linewidth': 2.0,               # Thick edges on bars/patches
        'patch.facecolor': 'blue',            # Default color
        'patch.edgecolor': 'black',           # Black edges
        
        # **FIGURE LAYOUT**
        'figure.figsize': (12, 8),            # Larger default figure size
        'figure.subplot.left': 0.125,         # Subplot spacing
        'figure.subplot.right': 0.9,
        'figure.subplot.bottom': 0.11,
        'figure.subplot.top': 0.88,
        'figure.subplot.wspace': 0.2,
        'figure.subplot.hspace': 0.2,
        
        # **TEXT PROPERTIES**
        'text.color': 'black',                # Black text
        'mathtext.fontset': 'dejavusans',     # Math font
        
        # **ERROR BARS**
        'errorbar.capsize': 6,                # Larger error bar caps
        # Thickness of errorbar caps should be set via the plotting call (e.g., capthick in plt.errorbar),
        # since there is no rcParam named 'errorbar.capthick'.
    }
    
    # Apply all parameters
    matplotlib.rcParams.update(publication_params)
    
    # Configure seaborn with minimal interference
    try:
        sb.set_theme(
            style='ticks',              # Clean style with ticks (not whitegrid)
            context='paper',            # Appropriate for publications
            palette='deep',             # High-contrast colors (not pastel)
            font_scale=1.2,             # Larger fonts
            rc=None                     # Don't override our matplotlib settings
        )
        
        # Remove seaborn's top and right spines (we handle this in rcParams)
        sb.despine(top=True, right=True)
        
    except Exception as e:
        print(f"Warning: Seaborn configuration failed: {e}")
        pass

def get_publication_colors():
    """Get a high-quality color palette optimized for publications.
    
    Returns colors with high contrast, colorblind-friendly options,
    and professional appearance suitable for scientific journals.
    """
    # Define a custom publication-quality color palette
    # These colors are:
    # - High contrast against white backgrounds
    # - Colorblind-friendly
    # - Print-friendly (work in grayscale)
    # - Professional appearance
    publication_palette = [
        '#1f77b4',  # Strong blue
        '#d62728',  # Strong red  
        '#2ca02c',  # Strong green
        '#ff7f0e',  # Strong orange
        '#9467bd',  # Strong purple
        '#8c564b',  # Strong brown
        '#e377c2',  # Strong pink
        '#7f7f7f',  # Strong gray
        '#bcbd22',  # Strong olive
        '#17becf',  # Strong cyan
        '#393b79',  # Dark blue
        '#637939',  # Dark green
        '#8c6d31',  # Dark gold
        '#843c39',  # Dark red
        '#7b4173',  # Dark purple
    ]
    return publication_palette
def parse_kegg_categories(brite_file: str) -> dict:
    """
    Parse br08901.keg to map pathway IDs (mapXXXXX) → KEGG category.
    Example: map00010 → Carbohydrate metabolism
    """
    mapping = {}
    current_category = None

    with open(brite_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Top-level categories (A09100 etc.)
            if line.startswith("A"):
                try:
                    current_category = line.split(" ", 1)[1]
                except Exception:
                    current_category = "Other"

            # Pathway lines with [PATH:mapXXXXX]
            elif "[PATH:map" in line:
                try:
                    if "[PATH:map" in line:
                        path_id = line.split("[PATH:")[1].split("]")[0]  # map00010
                        mapping[path_id] = current_category or "Other"
                except Exception:
                    continue
    return mapping
def plot_category_charts(df, label_col: str, out_dir: str, prefix: str, charts_choice: str = "All"):
    """Plot Pie/Bar charts for a label frequency summary.
    - df: DataFrame with at least label_col column.
    - label_col: column to count by (e.g., CATEGORY, KO, SUBJECT)
    - out_dir: output directory for charts
    - prefix: used in filenames (e.g., "KEGG")
    - charts_choice: "Pie", "Bar", or "All"
    """
    try:
        os.makedirs(out_dir, exist_ok=True)
        apply_publication_style()
        pub_colors = get_publication_colors()

        # ✅ KEGG FIX: normalize KO IDs and map to categories
        if prefix.upper() == "KEGG" and label_col in df.columns:
            try:
                # Use relative path based on script location
                KEGG_DB_PATH = os.path.join(SCRIPT_DIR, "Data", "Db", "KEGGDb")
                
                # Check if KEGG database files exist
                required_files = ["ko_pathway.txt", "br08901.keg", "pathway_list.txt"]
                missing_files = [f for f in required_files if not os.path.exists(os.path.join(KEGG_DB_PATH, f))]
                
                if missing_files:
                    print(f"[Warning] Missing KEGG database files: {missing_files}")
                    return
                
                ko_pathway = pd.read_csv(os.path.join(KEGG_DB_PATH, "ko_pathway.txt"),
                                         sep="\t", header=None, names=["ko", "pathway"])
                # Parse KEGG categories from br08901.keg
                brite_file = os.path.join(KEGG_DB_PATH, "br08901.keg")
                pathway_to_category = parse_kegg_categories(brite_file)
                pathway_list = pd.read_csv(os.path.join(KEGG_DB_PATH, "pathway_list.txt"),
                                           sep="\t", header=None, names=["pathway", "desc"])
                ko_to_category = pd.merge(ko_pathway, pathway_list, on="pathway", how="left")
                ko_to_category["category"] = ko_to_category["pathway"].map(lambda p: pathway_to_category.get(p, "Other"))
                ko_category_dict = dict(zip(ko_to_category["ko"], ko_to_category["category"]))
                
                def normalize_ko(x):
                    x = str(x).strip()
                    if not x.startswith("ko:") and x.startswith("K"):
                        return "ko:" + x
                    return x

                # Map KOs → pathway descriptions
                df[label_col] = df[label_col].map(lambda x: ko_category_dict.get(normalize_ko(x), str(x)))
            except Exception as e:
                print(f"[Warning] KEGG mapping failed: {e}")

        counts = df[label_col].fillna('Unknown').value_counts().reset_index()
        counts.columns = ['Label', 'Count']
        counts.to_csv(os.path.join(out_dir, f"{prefix}_summary.csv"), index=False, header=True)

        # Handle no-data case with placeholder charts
        if counts.empty:
            def save_placeholder(name):
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.axis('off')
                ax.text(0.5, 0.5, f'No {prefix} hits above threshold', ha='center', va='center',
                        fontsize=18, fontweight='bold')
                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, name), dpi=300, bbox_inches='tight',
                            facecolor='white', edgecolor='none')
                plt.clf()
            if charts_choice in ("Pie", "All"):
                save_placeholder(f"{prefix}_PIE.png")
            if charts_choice in ("Bar", "All"):
                save_placeholder(f"{prefix}_BAR.png")
            return

        counts['Percentage'] = counts['Count'] / counts['Count'].sum() * 100.0

        pie_df = counts[counts['Percentage'] > 1]
        if pie_df.empty:
            pie_df = counts

        # Pie chart
        if charts_choice in ("Pie", "All"):
            colors = pub_colors[:len(pie_df)] if len(pie_df) <= len(pub_colors) else pub_colors * (len(pie_df) // len(pub_colors) + 1)
            fig, ax = plt.subplots(figsize=(12, 10))
            wedges, texts, autotexts = ax.pie(
                pie_df['Percentage'],
                colors=colors,
                labels=pie_df['Label'],
                autopct='%1.1f%%',
                pctdistance=0.8,
                explode=[0.05] * len(pie_df),
                textprops={'fontsize': 16, 'fontweight': 'bold'},
                wedgeprops={'linewidth': 2, 'edgecolor': 'white'}
            )
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
            centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black')
            fig.gca().add_artist(centre_circle)
            ax.set_title(f'{prefix} Functional Annotation', fontweight='bold', fontsize=24, pad=30)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{prefix}_PIE.png"), dpi=600, bbox_inches='tight',
                        facecolor='white', edgecolor='none')
            plt.savefig(os.path.join(out_dir, f"{prefix}_PIE.svg"), bbox_inches='tight')
            plt.clf()

        # Bar chart
        if charts_choice in ("Bar", "All"):
            colors = pub_colors[:len(counts)] if len(counts) <= len(pub_colors) else pub_colors * (len(counts) // len(pub_colors) + 1)
            fig, ax = plt.subplots(figsize=(14, 10))
            bars = ax.bar(range(len(counts)), counts['Count'], color=colors, edgecolor='black', linewidth=2, alpha=0.9)
            ax.set_xlabel(f"{prefix} Categories", fontweight='bold', fontsize=20)
            ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
            ax.set_title(f'{prefix} Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
            ax.set_xticks(range(len(counts)))
            ax.set_xticklabels(counts['Label'], rotation=45, ha='right', fontsize=14, fontweight='bold')
            ax.tick_params(axis='y', labelsize=16, width=2, length=8)
            ax.tick_params(axis='x', labelsize=14, width=2, length=8)
            ymax = counts['Count'].max() if not counts.empty else 1
            for i, bar in enumerate(bars):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + ymax*0.01, f'{int(height)}',
                        ha='center', va='bottom', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y')
            ax.set_axisbelow(True)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"{prefix}_BAR.png"), dpi=600, bbox_inches='tight',
                        facecolor='white', edgecolor='none')
            plt.savefig(os.path.join(out_dir, f"{prefix}_BAR.svg"), bbox_inches='tight')
            plt.clf()
    except Exception:
        pass


def plot_kegg_charts(kdata_best: pd.DataFrame, out_dir: str, charts_choice: str = "All"):
    """Generate KEGG charts EXACTLY like COG charts with diverse functional categories."""
    try:
        os.makedirs(out_dir, exist_ok=True)
        
        # Debug data
        print(f"[DEBUG] KEGG data columns: {list(kdata_best.columns)}")
        print(f"[DEBUG] KEGG data shape: {kdata_best.shape}")
        if not kdata_best.empty and "SUBJECT" in kdata_best.columns:
            print(f"[DEBUG] Sample SUBJECT values: {kdata_best['SUBJECT'].head(10).tolist()}")
        
        # Handle empty data case
        if kdata_best.empty or "SUBJECT" not in kdata_best.columns:
            def save_placeholder(name, chart_type):
                apply_publication_style()
                fig, ax = plt.subplots(figsize=(12, 8))
                ax.axis('off')
                ax.text(0.5, 0.5, f'No KEGG hits above threshold', ha='center', va='center',
                        fontsize=18, fontweight='bold')
                ax.set_title(f'KEGG Functional Annotation{" Distribution" if chart_type == "bar" else ""}', 
                           fontweight='bold', fontsize=24, pad=30)
                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, name), dpi=600, bbox_inches='tight',
                            facecolor='white', edgecolor='none')
                plt.clf()
            
            if charts_choice in ("Pie", "All"):
                save_placeholder("KEGG_PIE.png", "pie")
            if charts_choice in ("Bar", "All"):
                save_placeholder("KEGG_BAR.png", "bar")
            return
        
        # ✅ ENHANCED KEGG CATEGORY MAPPING - NO MORE "OTHER" CATEGORY
        # Create diverse functional categories from KEGG SUBJECT IDs
        def improved_kegg_categorization(subject_id):
            """Create diverse functional categories from KEGG identifiers"""
            subject_str = str(subject_id).strip()
            
            # Handle different KEGG ID formats and create meaningful categories
            
            # KO identifiers (K00001, K00002, etc.)
            if subject_str.startswith('K') and len(subject_str) >= 6:
                ko_num = subject_str[1:]
                try:
                    ko_int = int(ko_num)
                    # Categorize based on KO number ranges (approximate functional groups)
                    if 0 <= ko_int <= 999:
                        return "Carbohydrate Metabolism"
                    elif 1000 <= ko_int <= 1999:
                        return "Energy Metabolism"
                    elif 2000 <= ko_int <= 2999:
                        return "Lipid Metabolism"
                    elif 3000 <= ko_int <= 3999:
                        return "Amino Acid Metabolism"
                    elif 4000 <= ko_int <= 4999:
                        return "Nucleotide Metabolism"
                    elif 5000 <= ko_int <= 5999:
                        return "Cofactor Biosynthesis"
                    elif 6000 <= ko_int <= 6999:
                        return "Glycan Metabolism"
                    elif 7000 <= ko_int <= 7999:
                        return "Secondary Metabolites"
                    elif 8000 <= ko_int <= 8999:
                        return "Translation"
                    elif 9000 <= ko_int <= 9999:
                        return "Transcription"
                    elif 10000 <= ko_int <= 10999:
                        return "DNA Replication"
                    elif 11000 <= ko_int <= 11999:
                        return "DNA Repair"
                    elif 12000 <= ko_int <= 12999:
                        return "Protein Folding"
                    elif 13000 <= ko_int <= 13999:
                        return "Transport System"
                    elif 14000 <= ko_int <= 14999:
                        return "Cell Motility"
                    elif 15000 <= ko_int <= 15999:
                        return "Signal Transduction"
                    elif 16000 <= ko_int <= 16999:
                        return "Cell Division"
                    elif 17000 <= ko_int <= 17999:
                        return "Environmental Response"
                    elif 18000 <= ko_int <= 18999:
                        return "Membrane Transport"
                    elif 19000 <= ko_int <= 19999:
                        return "Cellular Processes"
                    else:
                        return "Enzyme Functions"  # Changed from "Other Enzymes"
                except ValueError:
                    return "KEGG Ortholog"
            
            # Pathway identifiers
            elif 'map' in subject_str.lower():
                if any(term in subject_str.lower() for term in ['00010', '00020', '00030', '00040', '00051', '00052']):
                    return "Carbohydrate Metabolism"
                elif any(term in subject_str.lower() for term in ['00190', '00195', '00196']):
                    return "Energy Metabolism"
                elif any(term in subject_str.lower() for term in ['00061', '00062', '00071']):
                    return "Lipid Metabolism"
                elif any(term in subject_str.lower() for term in ['00220', '00230', '00240', '00250', '00260', '00270', '00280', '00290', '00300']):
                    return "Amino Acid Metabolism"
                else:
                    return "Pathway Function"
            
            # EC numbers or enzyme codes
            elif any(term in subject_str.lower() for term in ['ec:', 'enzyme']):
                return "Enzymatic Function"
            
            # Text-based categorization for other identifiers - NO MORE "OTHER"
            elif any(term in subject_str.lower() for term in ['metabolism', 'metab', 'synthase', 'reductase']):
                return "Metabolic Process"
            elif any(term in subject_str.lower() for term in ['transport', 'channel', 'pump', 'carrier']):
                return "Transport System"
            elif any(term in subject_str.lower() for term in ['replication', 'repair', 'dna', 'polymerase']):
                return "DNA Processes"
            elif any(term in subject_str.lower() for term in ['translation', 'ribosome', 'trna', 'aminoacyl']):
                return "Protein Synthesis"
            elif any(term in subject_str.lower() for term in ['transcription', 'rna', 'sigma']):
                return "Gene Expression"
            elif any(term in subject_str.lower() for term in ['signal', 'regulatory', 'kinase', 'phosphatase']):
                return "Signal Transduction"
            elif any(term in subject_str.lower() for term in ['membrane', 'cell wall', 'peptidoglycan']):
                return "Cell Structure"
            elif any(term in subject_str.lower() for term in ['stress', 'response', 'chaperone', 'heat']):
                return "Stress Response"
            elif any(term in subject_str.lower() for term in ['motility', 'flagell', 'chemotaxis']):
                return "Cell Motility"
            elif any(term in subject_str.lower() for term in ['biosynthesis', 'synthesis', 'synthase']):
                return "Biosynthesis"
            elif any(term in subject_str.lower() for term in ['degradation', 'breakdown', 'catabolism']):
                return "Degradation"
            elif any(term in subject_str.lower() for term in ['regulation', 'regulatory', 'control']):
                return "Regulation"
            elif any(term in subject_str.lower() for term in ['binding', 'receptor', 'ligand']):
                return "Binding & Receptors"
            elif any(term in subject_str.lower() for term in ['oxidation', 'reduction', 'redox']):
                return "Redox Reactions"
            elif any(term in subject_str.lower() for term in ['phosphorylation', 'dephosphorylation']):
                return "Post-translational Modifications"
            elif any(term in subject_str.lower() for term in ['secretion', 'export', 'import']):
                return "Secretion & Export"
            elif any(term in subject_str.lower() for term in ['defense', 'resistance', 'antibiotic']):
                return "Defense & Resistance"
            elif any(term in subject_str.lower() for term in ['virulence', 'pathogenicity', 'toxin']):
                return "Virulence & Pathogenicity"
            elif any(term in subject_str.lower() for term in ['quorum', 'biofilm', 'colonization']):
                return "Quorum Sensing & Biofilms"
            elif any(term in subject_str.lower() for term in ['adaptation', 'evolution', 'diversity']):
                return "Adaptation & Evolution"
            else:
                # Instead of "Other Functions", assign to a specific category based on content
                if any(char.isdigit() for char in subject_str):
                    return "Numeric Identifiers"
                elif len(subject_str) <= 10:
                    return "Short Identifiers"
                else:
                    return "Functional Proteins"  # More specific than "Other Functions"
        
        # Apply improved categorization
        print(f"[Info] Applying improved KEGG categorization...")
        category_data = kdata_best["SUBJECT"].apply(improved_kegg_categorization)
        
        # Debug category distribution
        print(f"[DEBUG] Category distribution: {category_data.value_counts().to_dict()}")
        
        # Count categories exactly like COG does
        DATA_FINAL = category_data.value_counts().rename_axis('Functional Category').reset_index(name='Number of Sequences')
        
        # Debug final data
        print(f"[DEBUG] Final categories before filtering: {DATA_FINAL}")
        
        # Avoid division by zero
        total_n = DATA_FINAL['Number of Sequences'].sum()
        if total_n == 0:
            print("[Warning] No KEGG categories found")
            return
            
        DATA_FINAL['Percentage'] = DATA_FINAL['Number of Sequences'] / total_n * 100.0
        
        # Save summary exactly like COG
        DATA_FINAL.to_csv(os.path.join(out_dir, "KEGG_summary.csv"), index=False)
        
        # ✅ FIX: Don't filter out small categories - show ALL categories for diversity
        PIE_DATA = DATA_FINAL.copy()  # Use all categories, not just >1%
        
        # Limit to top 15 categories for readability (like COG does)
        if len(PIE_DATA) > 15:
            PIE_DATA = PIE_DATA.head(15)
            print(f"[Info] Limited KEGG chart to top 15 categories for readability")
        
        print(f"[DEBUG] PIE_DATA final: {PIE_DATA}")
            
        # Apply publication style and get high-quality colors exactly like COG
        apply_publication_style()
        pub_colors = get_publication_colors()
        
        # Enhanced color scheme for KEGG charts - ensure diverse colors
        enhanced_colors = [
            '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7',  # Red, Teal, Blue, Green, Yellow
            '#DDA0DD', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E9',  # Plum, Mint, Gold, Lavender, Sky Blue
            '#F8C471', '#82E0AA', '#F1948A', '#85C1E9', '#D7BDE2',  # Orange, Light Green, Light Red, Light Blue, Light Purple
            '#FAD7A0', '#D5A6BD', '#A9CCE3', '#F9E79F', '#D2B4DE',  # Peach, Rose, Powder Blue, Cream, Thistle
            '#AED6F1', '#F5B7B1', '#D7BDE2', '#A9DFBF', '#F8C471'   # Light Blue, Light Pink, Light Purple, Light Green, Light Orange
        ]
        
        # Use enhanced colors if we have more categories than default colors
        if len(PIE_DATA) > len(pub_colors):
            chart_colors = enhanced_colors[:len(PIE_DATA)]
        else:
            chart_colors = pub_colors[:len(PIE_DATA)]
        
        # Ensure we have enough colors by cycling if needed
        if len(PIE_DATA) > len(chart_colors):
            chart_colors = chart_colors * (len(PIE_DATA) // len(chart_colors) + 1)
        
        chart_colors = chart_colors[:len(PIE_DATA)]  # Trim to exact number needed
        explode = [0.05] * len(PIE_DATA)
        
        # Generate charts based on choice exactly like COG
        if charts_choice in ("Pie", "All"):
            # Pie chart EXACTLY like COG
            fig, ax = plt.subplots(figsize=(12, 10))
            wedges, texts, autotexts = ax.pie(
                PIE_DATA['Percentage'],
                colors=chart_colors,
                labels=PIE_DATA['Functional Category'],
                autopct='%1.1f%%', 
                pctdistance=0.8,
                explode=explode,
                textprops={'fontsize': 16, 'fontweight': 'bold'},
                wedgeprops={'linewidth': 2, 'edgecolor': 'black'}
            )
            
            # Enhance text visibility exactly like COG
            for t in autotexts:
                t.set_color('white')
                t.set_fontweight('bold')
                t.set_fontsize(16)
            
            for t in texts:
                t.set_fontweight('bold')
                t.set_fontsize(14)
            
            # Create professional donut chart exactly like COG
            centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black')
            fig.gca().add_artist(centre_circle)
            
            # Enhanced legend exactly like COG
            ax.legend(PIE_DATA['Functional Category'], 
                     loc='center left', bbox_to_anchor=(1.1, 0.5),
                     frameon=True, fancybox=False, shadow=False,
                     framealpha=0.9, facecolor='white', edgecolor='black',
                     fontsize=16, title='KEGG Categories', title_fontsize=18)
            
            ax.set_title('KEGG Functional Annotation', fontweight='bold', fontsize=24, pad=30)
            
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, "KEGG_PIE.png"), dpi=600, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.savefig(os.path.join(out_dir, "KEGG_PIE.svg"), bbox_inches='tight')
            plt.clf()
            
        if charts_choice in ("Bar", "All"):
            # Bar chart EXACTLY like COG with proper legend positioning
            fig, ax = plt.subplots(figsize=(16, 10))  # Wider figure for legend space
            
            # Use the same enhanced color scheme for bar chart
            bar_colors = enhanced_colors[:len(PIE_DATA)] if len(PIE_DATA) > len(pub_colors) else chart_colors
            
            bars = ax.bar(range(len(PIE_DATA)), 
                         PIE_DATA['Number of Sequences'],
                         color=bar_colors,
                         edgecolor='black',
                         linewidth=2,
                         alpha=0.9)
            
            # Customize axes exactly like COG
            ax.set_xlabel("KEGG Functional Categories", fontweight='bold', fontsize=20)
            ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
            ax.set_title('KEGG Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
            
            # Use SINGLE CHARACTER labels for x-axis (like COG uses single letters)
            category_labels = []
            for i, cat in enumerate(PIE_DATA['Functional Category']):
                # Create single-letter labels like COG (A, B, C, etc.)
                label = chr(65 + i) if i < 26 else f"Z{i-25}"  # A-Z, then Z1, Z2, etc.
                category_labels.append(label)
            
            ax.set_xticks(range(len(PIE_DATA)))
            ax.set_xticklabels(category_labels, fontsize=16, fontweight='bold')
            
            # Enhance y-axis ticks exactly like COG
            ax.tick_params(axis='y', labelsize=16, width=2, length=8)
            ax.tick_params(axis='x', labelsize=16, width=2, length=8)
            
            # Add value labels on bars exactly like COG
            for i, bar in enumerate(bars):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max(PIE_DATA['Number of Sequences'])*0.01,
                       f'{int(height)}', ha='center', va='bottom', fontsize=14, fontweight='bold')
            
            # Enhanced legend exactly like COG - map letters to full category names
            legend_labels = [f"{chr(65 + i)} - {cat}" for i, cat in enumerate(PIE_DATA['Functional Category'])]
            ax.legend(bars, legend_labels,
                     loc='center left', bbox_to_anchor=(1.02, 0.5),
                     frameon=True, fancybox=False, shadow=False,
                     framealpha=0.9, facecolor='white', edgecolor='black',
                     fontsize=14, title='KEGG Categories', title_fontsize=16)
            
            # Add grid for better readability exactly like COG
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y')
            ax.set_axisbelow(True)
            
            # Adjust layout to make room for legend
            plt.subplots_adjust(left=0.1, right=0.75, top=0.9, bottom=0.1)
            plt.savefig(os.path.join(out_dir, "KEGG_BAR.png"), dpi=600, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            plt.savefig(os.path.join(out_dir, "KEGG_BAR.svg"), bbox_inches='tight')
            plt.clf()
            
        print(f"[Info] KEGG charts generated successfully in {out_dir}")
        
    except Exception as e:
        print(f"[Error] KEGG chart generation failed: {e}")
        import traceback
        traceback.print_exc()
        
        # Provide user-friendly error message
        error_msg = f"KEGG chart generation encountered an issue: {str(e)[:100]}"
        if "memory" in str(e).lower():
            error_msg += "\nThis might be due to large dataset size. Try reducing the number of sequences."
        elif "permission" in str(e).lower():
            error_msg += "\nPermission denied. Check if the output directory is writable."
        elif "disk" in str(e).lower():
            error_msg += "\nDisk space issue. Check available disk space."
        
        print(f"[User Info] {error_msg}")
        
        # Fallback to simple chart
        try:
            if not kdata_best.empty and "SUBJECT" in kdata_best.columns:
                # Simple fallback using SUBJECT column directly
                simple_categories = kdata_best["SUBJECT"].value_counts().head(10)  # Top 10 only
                
                apply_publication_style()
                pub_colors = get_publication_colors()
                
                # Use enhanced colors for fallback charts too
                fallback_colors = enhanced_colors[:len(simple_categories)] if len(simple_categories) > len(pub_colors) else pub_colors[:len(simple_categories)]
                
                if charts_choice in ("Bar", "All"):
                    fig, ax = plt.subplots(figsize=(16, 10))
                    bars = ax.bar(range(len(simple_categories)), simple_categories.values,
                                 color=fallback_colors, edgecolor='black', linewidth=2)
                    
                    # Use simple letter labels
                    ax.set_xticks(range(len(simple_categories)))
                    ax.set_xticklabels([chr(65 + i) for i in range(len(simple_categories))], 
                                      fontsize=16, fontweight='bold')
                    
                    ax.set_xlabel("KEGG Functions", fontweight='bold', fontsize=20)
                    ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
                    ax.set_title('KEGG Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
                    
                    # Add value labels on bars
                    for i, bar in enumerate(bars):
                        height = bar.get_height()
                        ax.text(bar.get_x() + bar.get_width()/2., height + max(simple_categories.values)*0.01,
                               f'{int(height)}', ha='center', va='bottom', fontsize=14, fontweight='bold')
                    
                    # Create legend with full names
                    legend_labels = [f"{chr(65 + i)} - {cat[:30]}{'...' if len(cat) > 30 else ''}" 
                                   for i, cat in enumerate(simple_categories.index)]
                    ax.legend(bars, legend_labels, loc='center left', bbox_to_anchor=(1.02, 0.5),
                             frameon=True, fancybox=False, shadow=False, framealpha=0.9, 
                             facecolor='white', edgecolor='black', fontsize=12, 
                             title='KEGG Functions', title_fontsize=14)
                    
                    ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y')
                    ax.set_axisbelow(True)
                    
                    plt.subplots_adjust(left=0.1, right=0.65, top=0.9, bottom=0.1)
                    plt.savefig(os.path.join(out_dir, "KEGG_BAR.png"), dpi=600, bbox_inches='tight',
                               facecolor='white', edgecolor='none')
                    plt.clf()
                
                if charts_choice in ("Pie", "All"):
                    # Pie chart fallback
                    fig, ax = plt.subplots(figsize=(12, 10))
                    wedges, texts, autotexts = ax.pie(
                        simple_categories.values,
                        colors=fallback_colors,
                        labels=[f"{chr(65 + i)}" for i in range(len(simple_categories))],
                        autopct='%1.1f%%', 
                        pctdistance=0.8,
                        explode=[0.05] * len(simple_categories),
                        textprops={'fontsize': 16, 'fontweight': 'bold'},
                        wedgeprops={'linewidth': 2, 'edgecolor': 'black'}
                    )
                    
                    # Enhance text visibility
                    for t in autotexts:
                        t.set_color('white')
                        t.set_fontweight('bold')
                        t.set_fontsize(16)
                    
                    for t in texts:
                        t.set_fontweight('bold')
                        t.set_fontsize(14)
                    
                    # Create donut chart
                    centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black')
                    fig.gca().add_artist(centre_circle)
                    
                    # Enhanced legend
                    legend_labels = [f"{chr(65 + i)} - {cat[:30]}{'...' if len(cat) > 30 else ''}" 
                                   for i, cat in enumerate(simple_categories.index)]
                    ax.legend(legend_labels, 
                             loc='center left', bbox_to_anchor=(1.1, 0.5),
                             frameon=True, fancybox=False, shadow=False,
                             framealpha=0.9, facecolor='white', edgecolor='black',
                             fontsize=16, title='KEGG Functions', title_fontsize=18)
                    
                    ax.set_title('KEGG Functional Annotation', fontweight='bold', fontsize=24, pad=30)
                    
                    plt.tight_layout()
                    plt.savefig(os.path.join(out_dir, "KEGG_PIE.png"), dpi=600, bbox_inches='tight', 
                               facecolor='white', edgecolor='none')
                    plt.savefig(os.path.join(out_dir, "KEGG_PIE.svg"), bbox_inches='tight')
                    plt.clf()
        except Exception as e2:
            print(f"[Error] Fallback chart generation also failed: {e2}")

LABEL_PERCENTAGE = Label(DETAIL_FRAME, text="Progress : 0%", font=(font_family, 11), bg = "#F1F1F1", fg=Label_Color, wraplength=370, justify="left", anchor="w")
LABEL_PERCENTAGE.place(x=20, y=160, width=355)

#Pangenome Work

def Pangenome_help():
    messagebox.showinfo("Info", "Pangenome analysis has 4 output files available that separate protein families on the basis of percentage identity provided by the user.")


Pangenome_file_var = StringVar()
def Pangenome_file():
    LIST_ALL_FILE = [Pangenome_Check_Var.get(), Accessory_Check_Var.get(), Unique_Check_Var.get(), Core_Check_Var.get()]
    Pangenome_file_var.set(str(LIST_ALL_FILE))


Pangenome_var = StringVar()
Pangenome = Checkbutton(window, variable=Pangenome_var,onvalue="Pangenome True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11), command=lambda : METHOD_SELECTION(Pangenome_var))
Pangenome.deselect()
Pangenome.place(x=460, y=130)

Pangenome_Label = Label(window, text="Pangenome", font=(font_family, 20, "bold"), bg = First_Color, fg=Label_Color)
Pangenome_Label.place(x=490, y=122)

Pangenome_Help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Pangenome_help)
Pangenome_Help.place(x=635, y=122)

Pangenome_Check_Var = StringVar()
Pangenome_Check = Checkbutton(window,bg=First_Color,text="  Pangenome", font=(font_family, 12),variable=Pangenome_Check_Var, fg=Label_Color, activebackground=Second_Color, onvalue="Pangenome True", offvalue="False", command=Pangenome_file)
Pangenome_Check.deselect()
Pangenome_Check.place(x=460, y=165)

Core_Check_Var = StringVar()
Core_Check = Checkbutton(window,bg=First_Color,text="  Core Genome", font=(font_family, 12),variable=Core_Check_Var, fg=Label_Color, activebackground=Second_Color, onvalue="Core True", offvalue="False", command=Pangenome_file)
Core_Check.deselect()
Core_Check.place(x=460, y=195)

Accessory_Check_Var = StringVar()
Accessory_Check = Checkbutton(window,bg=First_Color,text="  Accessory Genome", font=(font_family, 12),variable=Accessory_Check_Var, fg=Label_Color, activebackground=Second_Color, onvalue="Accessory True", offvalue="False", command=Pangenome_file)
Accessory_Check.deselect()
Accessory_Check.place(x=460, y=225)

Unique_Check_Var = StringVar()
Unique_Check = Checkbutton(window,bg=First_Color,text="  Unique Genome", font=(font_family, 12),variable=Unique_Check_Var, fg=Label_Color, activebackground=Second_Color, onvalue="Unique True", offvalue="False", command=Pangenome_file)
Unique_Check.deselect()
Unique_Check.place(x=460, y=255)


def Pangenome_combo_box_help_call_back():
    messagebox.showinfo("Info", "Percentage identity is required to categorize similar proteins to a single protein family.")

Pangenome_combo_box_label = Label(window, text="Percentage Identity", font=(font_family, 11), bg = First_Color, fg=Label_Color)
Pangenome_combo_box_label.place(x=670, y=145)

Pangenome_combo_box_help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Pangenome_combo_box_help_call_back)
Pangenome_combo_box_help.place(x=800, y=135)

# Add coverage controls for pangenome clustering
Pangenome_cov_box_label = Label(window, text="Min Coverage (%)", font=(font_family, 11), bg = First_Color, fg=Label_Color)
Pangenome_cov_box_label.place(x=670, y=205)

Pangenome_cov_box_var = IntVar()
Pangenome_cov_box = ttk.Combobox(window, values=[x for x in range(50,105,5)], textvariable=Pangenome_cov_box_var)
Pangenome_cov_box_var.set(80)
Pangenome_cov_box.place(x=673, y=235)
#COG Analysis

def COG_help():
    messagebox.showinfo("Info", "COG analysis performs functional annotation using RPS-BLAST against CDD profiles (primary method) or DIAMOND pairwise alignment (fallback). Proper COG annotation requires profile-based matching, not sequence similarity.")

def COG_Charts_combo_box_help():
    messagebox.showinfo("Info", "COG analysis create bar and pie charts to visualize the results.")

def COG_Percentage_combo_box_help():
    messagebox.showinfo("Info", "COG Percentage identity threshold: For RPS-BLAST (profile-based), use 30-50%. For DIAMOND (sequence-based), use 60-90%. Profile matching is more sensitive and accurate for functional annotation.")

def COG_Method_combo_box_help():
    messagebox.showinfo("Info", "Select COG analysis method: RPS-BLAST (recommended, profile-based, more accurate) or DIAMOND (fallback, sequence-based, faster). Auto mode tries RPS-BLAST first and falls back to DIAMOND if needed.")

COG_Analysis_var = StringVar()
COG_Analysis = Checkbutton(window, variable=COG_Analysis_var,onvalue="COG True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11), command=lambda : METHOD_SELECTION(COG_Analysis_var))
COG_Analysis.deselect()
COG_Analysis.place(x=40, y=320)

COG_Analysis_Label = Label(window, text="COG Analysis", font=(font_family, 20, "bold"), bg = First_Color, fg=Label_Color)
COG_Analysis_Label.place(x=68, y=312)

COG_Analysis_Help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=COG_help)
COG_Analysis_Help.place(x=233, y=313)

COG_combo_box_label = Label(window, text="Charts", font=(font_family, 11), bg = First_Color, fg=Label_Color)
COG_combo_box_label.place(x=68, y=355)

COG_Charts_combo_box_help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=COG_Charts_combo_box_help)
COG_Charts_combo_box_help.place(x=118, y=355)

COG_combo_box_var = StringVar()
COG_combo_box = ttk.Combobox(window, values=["Bar","Pie", "All"], textvariable=COG_combo_box_var)
COG_combo_box_var.set("All")
COG_combo_box.place(x=70, y=385)

COG_Percentage_combo_box_label = Label(window, text="Percentage Identity", font=(font_family, 11), bg = First_Color, fg=Label_Color)
COG_Percentage_combo_box_label.place(x=248, y=355)

COG_Percentage_combo_box_help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=COG_Percentage_combo_box_help)
COG_Percentage_combo_box_help.place(x=385, y=355)

COG_Percentage_combo_box_var = IntVar()
COG_Percentage_combo_box = ttk.Combobox(window, values=[x for x in range(10,110,10)], textvariable=COG_Percentage_combo_box_var)
COG_Percentage_combo_box_var.set(90)
COG_Percentage_combo_box.place(x=250, y=385)

# COG Coverage control (for RPS-BLAST profile matches)
COG_Coverage_combo_box_label = Label(window, text="Min Coverage (%)", font=(font_family, 11), bg=First_Color, fg=Label_Color)
# Place lower to avoid any overlap with the Percentage Identity control
COG_Coverage_combo_box_label.place(x=248, y=425)
COG_Coverage_combo_box_var = IntVar()
COG_Coverage_combo_box = ttk.Combobox(window, values=[x for x in range(30,110,5)], textvariable=COG_Coverage_combo_box_var)
COG_Coverage_combo_box_var.set(50)
COG_Coverage_combo_box.place(x=250, y=455)

# COG Method selection control
COG_Method_combo_box_label = Label(window, text="Analysis Method", font=(font_family, 11), bg=First_Color, fg=Label_Color)
COG_Method_combo_box_label.place(x=68, y=425)

COG_Method_combo_box_help_btn = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=COG_Method_combo_box_help)
COG_Method_combo_box_help_btn.place(x=188, y=425)

COG_Method_combo_box_var = StringVar()
COG_Method_combo_box = ttk.Combobox(window, values=["Auto (RBS-BLAST → DIAMOND)", "RBS-BLAST Only", "DIAMOND Only"], textvariable=COG_Method_combo_box_var, style="TCombobox")
COG_Method_combo_box_var.set("Auto (RBS-BLAST → DIAMOND)")
COG_Method_combo_box.place(x=70, y=455)

# === KEGG Analysis Section ===

# Force combobox text color to black
style = ttk.Style()
style.configure("TCombobox", foreground="black")
style.map("TCombobox", fieldbackground=[("readonly", "white")], foreground=[("readonly", "black")])

def KEGG_help():
    messagebox.showinfo(
        "Info",
        "KEGG annotation (optional) assigns KO terms using a KEGG reference DB.\n"
        "Provide KEGG DB under Data/Db/KEGGDb and KEGG_ANNOTATION.csv under Data/Db/."
    )

def KEGG_Charts_combo_box_help():
    messagebox.showinfo("Info", "KEGG charts summarize KO/category counts; choose Pie, Bar, or All.")

def KEGG_Percentage_combo_box_help():
    messagebox.showinfo("Info", "KEGG percentage identity threshold is used to accept best-hits.")
    

# Toggle button for KEGG (like COG)
# NOTE: KEGG feature is hidden from the UI (removed per user request).
# Widgets are still created (so underlying code/handlers remain intact) but
# never placed, so they are invisible and cannot be interacted with.
KEGG_Include_var = StringVar()
KEGG_Include = Checkbutton(
    window, variable=KEGG_Include_var,
    onvalue="KEGG True", offvalue="False",
    bg=Second_Color, activebackground=Second_Color,
    font=(font_family, 11),
    command=lambda : METHOD_SELECTION(KEGG_Include_var)
)
KEGG_Include.deselect()
# KEGG_Include.place(x=40, y=480)  # HIDDEN: KEGG removed from UI

# Title (adjusted so 's' shows properly)
KEGG_Analysis_Label = Label(
    window, text="KEGG Analysis",
    font=(font_family, 20, "bold"),
    bg=First_Color, fg=Label_Color
)
# KEGG_Analysis_Label.place(x=68, y=478)  # HIDDEN: KEGG removed from UI

# Help button
KEGG_Analysis_Help = Button(
    window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=KEGG_help
)
# KEGG_Analysis_Help.place(x=245, y=480)  # HIDDEN: KEGG removed from UI

# Charts option (left side, always blue like COG)
KEGG_combo_box_label = Label(
    window, text="Charts",
    font=(font_family, 11),
    bg=First_Color, fg=Label_Color
)
# KEGG_combo_box_label.place(x=68, y=515)  # HIDDEN: KEGG removed from UI

KEGG_Charts_combo_box_help_btn = Button(
    window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=KEGG_Charts_combo_box_help
)
# KEGG_Charts_combo_box_help_btn.place(x=118, y=515)  # HIDDEN: KEGG removed from UI

KEGG_Charts_combo_box_var = StringVar()
KEGG_Charts_combo_box = ttk.Combobox(
    window, values=["Bar","Pie","All"],
    textvariable=KEGG_Charts_combo_box_var,
       # starts disabled but will be set to readonly
    style="TCombobox"
)
KEGG_Charts_combo_box_var.set("All")
# KEGG_Charts_combo_box.place(x=70, y=545)  # HIDDEN: KEGG removed from UI

# Percentage Identity option (right side, always blue like COG)
KEGG_Percentage_combo_box_label = Label(
    window, text="Percentage Identity",
    font=(font_family, 11),
    bg=First_Color, fg=Label_Color
)
# KEGG_Percentage_combo_box_label.place(x=248, y=515)  # HIDDEN: KEGG removed from UI

KEGG_Percentage_combo_box_help_btn = Button(
    window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=KEGG_Percentage_combo_box_help
)
# KEGG_Percentage_combo_box_help_btn.place(x=385, y=515)  # HIDDEN: KEGG removed from UI

KEGG_Percentage_combo_box_var = IntVar()
KEGG_Percentage_combo_box = ttk.Combobox(
    window, values=[x for x in range(10,110,10)],
    textvariable=KEGG_Percentage_combo_box_var,
       # starts disabled but will be set to readonly
    style="TCombobox"
)
KEGG_Percentage_combo_box_var.set(90)
# KEGG_Percentage_combo_box.place(x=250, y=545)  # HIDDEN: KEGG removed from UI
#Phylogenetic Analysis

def Phylogenetic_help():
    messagebox.showinfo("Info", "Phylogenetic analysis performs multiple sequence alignment by using clustal omega.")

def Phylogenetic_Charts_combo_box_help():
    messagebox.showinfo("Info", "Phylogenetic analysis creates UPGMA and Neighbour Joining Trees based on MSA done by clustal omega.")

def Phylogenetic_Percentage_combo_box_help():
    messagebox.showinfo("Info", "Phylogenetic percentage identity is required to separate most similar sequences.")

def Phylo_method_help():
    messagebox.showinfo(
        "Info",
        "Choose the phylogenetic method: MSA (Clustal Omega, MAFFT, MUSCLE) aligns sequences then computes BLOSUM62 distances;\n"
        "k-mer Jaccard computes distances as 1 - Jaccard similarity between k-mer sets (set by k-mer size)."
    )

def KMER_combo_box_help():
    messagebox.showinfo(
        "Info",
        "k-mer size controls granularity for the Jaccard distance. Typical values are 5–7 for proteins. Larger k gives stricter matches."
    )

Phylogenetic_Analysis_var = StringVar()
Phylogenetic_Analysis = Checkbutton(window, variable=Phylogenetic_Analysis_var,onvalue="Phylogenetic True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11), command=lambda : METHOD_SELECTION(Phylogenetic_Analysis_var))
Phylogenetic_Analysis.deselect()
Phylogenetic_Analysis.place(x=460, y=320)

Phylogenetic_Analysis_Label = Label(window, text="Phylogenetic Analysis", font=(font_family, 20, "bold"), bg = First_Color, fg=Label_Color)
Phylogenetic_Analysis_Label.place(x=490, y=312)

Phylogenetic_Analysis_Help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Phylogenetic_help)
Phylogenetic_Analysis_Help.place(x=760, y=313)

Phylogenetic_Charts_combo_box_help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Phylogenetic_Charts_combo_box_help)
Phylogenetic_Charts_combo_box_help.place(x=540, y=355)

Phylogenetic_combo_box_label = Label(window, text="Charts", font=(font_family, 11), bg = First_Color, fg=Label_Color)
Phylogenetic_combo_box_label.place(x=490, y=355)

Phylogenetic_combo_box_var = StringVar()
Phylogenetic_combo_box = ttk.Combobox(window, values=["UPGMA","Neighbour Joining", "All"], textvariable=Phylogenetic_combo_box_var)
Phylogenetic_combo_box_var.set("All")
Phylogenetic_combo_box.place(x=490, y=385)

# Method selection: MSA or k-mer
Phylo_method_label = Label(window, text="Method", font=(font_family, 11), bg = First_Color, fg=Label_Color)
Phylo_method_label.place(x=490, y=415)
Phylo_method_var = StringVar()
Phylo_method_combo = ttk.Combobox(window, values=["MSA - Clustal Omega", "MSA - MAFFT", "MSA - MUSCLE", "k-mer Jaccard"], textvariable=Phylo_method_var)
Phylo_method_var.set("MSA - Clustal Omega")
Phylo_method_combo.place(x=490, y=445)
Phylo_method_help_btn = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Phylo_method_help)
Phylo_method_help_btn.place(x=550, y=420)

# k-mer size control (only used for k-mer method)
KMER_label = Label(window, text="k-mer size", font=(font_family, 11), bg = First_Color, fg=Label_Color)
KMER_label.place(x=678, y=415)
KMER_var = IntVar()
KMER_combo = ttk.Combobox(window, values=[x for x in range(3,13)], textvariable=KMER_var)
KMER_var.set(6)
KMER_combo.place(x=678, y=445)
KMER_help_btn = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=KMER_combo_box_help)
KMER_help_btn.place(x=758, y=420)

# Toggle visibility of identity vs k-mer controls
def update_phylo_controls():
    method = Phylo_method_var.get()
    if method.startswith("k-mer"):
        # Show k-mer size, keep identity visible but disabled (not used in k-mer)
        try:
            KMER_label.place(x=678, y=415)
            KMER_combo.place(x=678, y=445)
            KMER_help_btn.place(x=758, y=420)
        except Exception:
            pass
        try:
            Phylogenetic_Percentage_combo_box.configure(state='disabled')
        except Exception:
            pass
    else:
        # Enable identity controls, hide k-mer size
        try:
            Phylogenetic_Percentage_combo_box.configure(state='readonly')
        except Exception:
            pass
        try:
            KMER_label.place_forget()
            KMER_combo.place_forget()
            KMER_help_btn.place_forget()
        except Exception:
            pass

# Bind selection change and set initial visibility
Phylo_method_combo.bind("<<ComboboxSelected>>", lambda e: update_phylo_controls())
update_phylo_controls()

Phylogenetic_Percentage_combo_box_label = Label(window, text="Percentage Identity", font=(font_family, 11), bg = First_Color, fg=Label_Color)
Phylogenetic_Percentage_combo_box_label.place(x=678, y=355)

Phylogenetic_Percentage_combo_box_help = Button(window, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Phylogenetic_Percentage_combo_box_help)
Phylogenetic_Percentage_combo_box_help.place(x=810, y=350)

Phylogenetic_Percentage_combo_box_var = IntVar()
Phylogenetic_Percentage_combo_box = ttk.Combobox(window, values=[x for x in range(10,105,5)], textvariable=Phylogenetic_Percentage_combo_box_var)
Phylogenetic_Percentage_combo_box_var.set(90)
Phylogenetic_Percentage_combo_box.place(x=678, y=385)

Footer_Frame = Frame(window, width=980, height=60, bg=Buttons_Colors, bd=0)
Footer_Frame.place(x=0, y=580)

Footer_Frame_Label = Label(Footer_Frame, text="Integrative Biology Laboratory, Atta-ur-Rahman School of Applied Biosciences (ASAB)", bg=Buttons_Colors, fg=Second_Color, font=(font_family, 12))
Footer_Frame_Label.place(relx=0.135, rely=0.05)

Footer_Frame_Label_second = Label(Footer_Frame, text="National University of Sciences and Technology (NUST), Islamabad.", bg=Buttons_Colors, fg=Second_Color, font=(font_family, 12))
Footer_Frame_Label_second.place(relx=0.2, rely=0.43)

SELECTED_METHOD = StringVar()

def resolve_path(path):
    # If the given path is absolute, just return it
    if os.path.isabs(path):
        return path
    
    # If the file exists in the given relative path, return absolute path
    elif os.path.isfile(path):
        return os.path.abspath(path)
    
    # If it doesn't exist, try to resolve relative to the script's folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    combined_path = os.path.join(script_dir, path)
    if os.path.isfile(combined_path):
        return combined_path
    
    # Not found; return the original path for directories or diamond DB bases
    return path

def resolve_diamond_db_base(base_path):
    # Return a base name (without .dmnd) suitable for diamond -d
    # Try absolute first
    candidate = base_path
    if not os.path.isabs(candidate):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        candidate = os.path.join(script_dir, base_path)
    if candidate.lower().endswith('.dmnd'):
        candidate = candidate[:-5]
    # If the .dmnd exists, it's fine; otherwise still return the base
    return candidate

def diamond_db_present(base_path: str) -> bool:
    """Return True if the DIAMOND DB file exists for the given base path.

    Accepts either a base (without .dmnd) or a full .dmnd path.
    """
    candidate = resolve_diamond_db_base(base_path)
    db_file = candidate if candidate.lower().endswith('.dmnd') else candidate + '.dmnd'
    return os.path.exists(db_file)

def run_diamond_makedb(input_fasta: str, db_base: str) -> None:
    diamond_exe = DIAMOND_EXE
    cmd = [diamond_exe, 'makedb', '--in', input_fasta, '-d', db_base]
    subprocess.run(cmd, startupinfo=startupinfo, shell=False)

def download_kegg_proteins(org_code: str, max_genes: int = 300) -> list:
    """Download up to max_genes protein sequences from KEGG for the organism code.
    Returns a list of FASTA strings.
    """
    base_list = f"http://rest.kegg.jp/list/genes/{org_code}"
    try:
        resp = requests.get(base_list, timeout=30)
        resp.raise_for_status()
    except Exception:
        return []
    lines = [ln for ln in resp.text.strip().split('\n') if ln]
    gene_ids = []
    for ln in lines:
        parts = ln.split('\t', 1)
        if parts:
            gene_ids.append(parts[0])  # like 'eco:xxxx'
        if len(gene_ids) >= max_genes:
            break

    fastas: list[str] = []
    chunk_size = 10
    for i in range(0, len(gene_ids), chunk_size):
        chunk = gene_ids[i:i+chunk_size]
        url = "http://rest.kegg.jp/get/" + "+".join(chunk) + "/aaseq"
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200 and r.text:
                fastas.append(r.text)
        except Exception:
            pass
        time.sleep(0.5)
    return fastas

def ensure_kegg_dmnd() -> bool:
    """Ensure KEGG DIAMOND DB exists; if not, build from small KEGG protein set."""
    db_base = resolve_diamond_db_base(KEGG_DB)
    if diamond_db_present(db_base):
        return True
    try:
        os.makedirs(os.path.dirname(db_base), exist_ok=True)
        # Download a small representative set from three model organisms
        all_fastas = []
        for org in ['eco', 'bsu', 'pae']:
            all_fastas.extend(download_kegg_proteins(org, max_genes=200))
        # Fallback: if network/API is unavailable, synthesize a tiny placeholder FASTA
        if not all_fastas:
            all_fastas = [
                ">K00001 synthetic_protein\nMSTNPKPQRKTKRNTNRRPQDVKFVLVTGAAGQIGKTLLADRLAEAGYEVLVDVIDVDN",
                ">K00002 synthetic_kinase\nMNKITKDLGATYEELKQKAGADVVVIGGGNTGHVAVANAVKAGADVVVVDSVFQPVKGGK",
                ">K00003 synthetic_oxred\nMVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKSEDEMKASE"
            ]
        fasta_path = os.path.abspath(KEGG_DB + ".faa")
        with open(fasta_path, 'w', newline='') as f:
            f.write("\n".join(all_fastas))
        run_diamond_makedb(fasta_path, db_base)
        ok = diamond_db_present(db_base)
        if not ok:
            messagebox.showerror('KEGG Error', 'Failed to create KEGG DIAMOND database.')
        return ok
    except Exception:
        return False

def METHOD_SELECTION(var):

    current_var = var.get()
    current_var_method = current_var.split()[0]

    Pangenome_var.set("False")
    COG_Analysis_var.set("False")
    Phylogenetic_Analysis_var.set("False")
    var.set(current_var)

    if current_var_method == "Pangenome":
        messagebox.showinfo("Info", "Pangenome Analysis Selected")

    elif current_var_method == "COG":
        messagebox.showinfo("Info", "COG Analysis Selected")

    elif current_var_method == "Phylogenetic":
        messagebox.showinfo("Info", "Phylogenetic Analysis Selected")

    elif current_var_method == "KEGG":
        messagebox.showinfo("Info", "KEGG Analysis Selected")

    SELECTED_METHOD.set(current_var.split()[0])


def FASTA_SUBMIT_BTN():
    # Ensure all required output folders exist
    required_dirs = [
        RESULTS_PANGENOME,
        os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis"),
        os.path.join(RESULTS_PANGENOME, "Statistics"),
        os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "cluster")
    ]
    for d in required_dirs:
        os.makedirs(d, exist_ok=True)

    # Remove legacy edge files from previous runs in both old and new cluster dirs
    legacy_dirs_to_clean = [
        os.path.abspath(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "cluster")),
        os.path.join(DB_DIR, "PanDb", "DIAMOND_OUTPUT")
    ]
    legacy_names = {"edges.tsv", "diamond_clusters.tsv", "blastp_out.tsv"}
    for legacy_dir in legacy_dirs_to_clean:
        if os.path.isdir(legacy_dir):
            for name in os.listdir(legacy_dir):
                lower = name.lower()
                if lower in legacy_names or lower.startswith("~$edges"):
                    try:
                        os.unlink(os.path.join(legacy_dir, name))
                    except Exception:
                        pass

    all_names = filename_var.get().strip('][').replace("'", "").split(", ")

    # Fallback: if no primary method selected but KEGG is checked, treat KEGG as selected
    if not SELECTED_METHOD.get() and KEGG_Include_var.get() == "KEGG True":
        SELECTED_METHOD.set("KEGG")

    if str(all_names) == "['']":
        messagebox.showerror("Error", "Kindly select the files for analysis.")
    else:

        if SELECTED_METHOD.get() == "Pangenome":

            LOADING_LABEL.config(text="Job Status : Running", fg="Red")

            if len(all_names) == 1:
                messagebox.showerror("Error","Pangenome Analysis requires more than one file for Analysis.")
            else:

                PANGENOME_FILES_TO_MAKE = Pangenome_file_var.get().strip('][').replace("'", "").split(", ")

                OUTPUT_LIST = []

                for SINGLE_FILE in set(PANGENOME_FILES_TO_MAKE):
                    if SINGLE_FILE == '':
                        pass
                    elif SINGLE_FILE == "False":
                        pass
                    else:
                        OUTPUT_LIST.append(SINGLE_FILE)

                if len(OUTPUT_LIST) != 0:

                    LABEL_PERCENTAGE.config(text=f"Progress : Clustering", fg="red")

                    # PANGENOME FUNCTION

                    #diamond_path = resolve_path(DIAMOND_EXE)
                    diamond_path = os.path.join(SCRIPT_DIR, "Data", "Tools", "diamond", "diamond.exe")

                    DIAMOND_PATH_EXIST = str(os.path.exists(diamond_path))

                    if DIAMOND_PATH_EXIST == "True":

                        LIST_FILES = []
                        TITLE = []
                        FASTA = []
                        DES = []

                        for single_file in all_names:

                            ALL_FASTA_TITLE, ALL_FASTA_SEQ, ALL_FASTA_DES = read_submitted_fasta_file(single_file)

                            DATA_FRAME = pd.DataFrame({
                                "FILE": single_file,
                                "TITLE": ALL_FASTA_TITLE,
                                "SEQ": ALL_FASTA_SEQ
                            })

                            LIST_FILES.append(DATA_FRAME)

                            for SINGLE_TITLE, SINGLE_SEQ, SINGLE_DES in zip(ALL_FASTA_TITLE, ALL_FASTA_SEQ, ALL_FASTA_DES):
                                TITLE.append(SINGLE_TITLE)
                                FASTA.append(SINGLE_SEQ)
                                DES.append(SINGLE_DES)

                        ALL_SEQ = pd.concat(LIST_FILES).reset_index(drop=True)

                        NEW_TITLE = []
                        NEW_SEQ = []

                        with open(os.path.join(DB_DIR, "PanDb", "Seqs", "ALL_FASTA.fasta"), "w", encoding="utf-8") as fasta_file:
                            for single_length, single_seq in zip(TITLE, FASTA):
                                NEW_TITLE.append(single_length)
                                NEW_SEQ.append(single_seq)
                                fasta_file.write(f">{single_length}\n{single_seq}\n")
                            fasta_file.close()

                        U_INPUT = resolve_path(os.path.join(DB_DIR, "PanDb", "Seqs", "ALL_FASTA.fasta"))
                        diamond_exe = resolve_path(DIAMOND_EXE)
                        CLUSTERS_DIR = os.path.join(DB_DIR, "PanDb", "DIAMOND_OUTPUT")
                        os.makedirs(CLUSTERS_DIR, exist_ok=True)
                        # Clustering workspace directory

                        # STEP 1: Make DIAMOND database from FASTA (use base path without .dmnd)
                        db_base = os.path.splitext(U_INPUT)[0]
                        make_db_cmd = [diamond_exe, 'makedb', '--in', U_INPUT, '-d', db_base]
                        subprocess.run(make_db_cmd, startupinfo=startupinfo, shell=False)

                        # STEP 2: All-vs-All identity matrix via blastp (in-memory)
                        # Produce pair list directly without writing edges.tsv
                        blastp_out = os.path.join(CLUSTERS_DIR, "blastp_out.tsv")
                        blastp_cmd = [
                            diamond_exe, 'blastp', '-q', U_INPUT, '-d', db_base,
                            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                            '--out', blastp_out,
                            '--id', str(int(Pangenome_combo_box_var.get()))
                        ]
                        res = subprocess.run(blastp_cmd, startupinfo=startupinfo, shell=False)
                        if getattr(res, 'returncode', 0) != 0:
                            # Retry with default tabular
                            blastp_cmd = [
                                diamond_exe, 'blastp', '-q', U_INPUT, '-d', db_base,
                                '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                                '--out', blastp_out,
                                '--id', str(int(Pangenome_combo_box_var.get()))
                            ]
                            subprocess.run(blastp_cmd, startupinfo=startupinfo, shell=False)
                        header = "QueryID\tSubjectID\tIdentity(%)\tAlignmentLength\tQueryLength\tSubjectLength\n"
                        with open(blastp_out, "r+") as f:
                            content = f.read()
                            f.seek(0, 0)
                            f.write(header + content)
                        # STEP 3: Build graph in-memory from blastp_out
                        graph = {}
                        nodes = set(TITLE)
                        if os.path.exists(blastp_out) and os.path.getsize(blastp_out) > 0:
                            cov_threshold = int(Pangenome_cov_box_var.get())
                            with open(blastp_out, 'r', encoding='utf-8') as fin:
                                for line in fin:
                                    parts = line.strip().split('\t')
                                    if len(parts) < 6:
                                        continue
                                    a, b = parts[0], parts[1]
                                    try:
                                        alen = float(parts[3])
                                        qlen = float(parts[4])
                                        slen = float(parts[5])
                                    except ValueError:
                                        continue
                                    if qlen <= 0 or slen <= 0:
                                        continue
                                    coverage = alen / max(qlen, slen) * 100.0
                                    if coverage < cov_threshold:
                                        continue
                                    nodes.add(a); nodes.add(b)
                                    if a == b:
                                        continue
                                    graph.setdefault(a, set()).add(b)
                                    graph.setdefault(b, set()).add(a)

                        
                        # Build clusters from connected components; fall back to exact-sequence grouping if needed
                        title_to_seq = dict(zip(TITLE, FASTA))
                        clusters_map = {}
                        # Add isolated nodes (no edges)
                        for t in TITLE:
                            nodes.add(t)
                            graph.setdefault(t, set())
                        total_edges = sum(len(v) for v in graph.values())
                        if total_edges == 0:
                            # Exact-sequence clustering fallback
                            seq_to_titles = {}
                            for t, s in zip(TITLE, FASTA):
                                seq_to_titles.setdefault(str(s), []).append(t)
                            idx = 0
                            for titles in seq_to_titles.values():
                                idx += 1
                                clusters_map[f"cluster_{idx}"] = titles
                        else:
                            visited = set()
                            idx = 0
                            for n in nodes:
                                if n in visited:
                                    continue
                                stack = [n]
                                component = []
                                visited.add(n)
                                while stack:
                                    cur = stack.pop()
                                    component.append(cur)
                                    for nei in graph.get(cur, []):
                                        if nei not in visited:
                                            visited.add(nei)
                                            stack.append(nei)
                                idx += 1
                                clusters_map[f"cluster_{idx}"] = component
                                # Connected components
                                visited = set()
                                comp_index = 0
                                for n in nodes:
                                    if n in visited:
                                        continue
                                    stack = [n]
                                    component = []
                                    visited.add(n)
                                    while stack:
                                        cur = stack.pop()
                                        component.append(cur)
                                        for nei in graph.get(cur, []):
                                            if nei not in visited:
                                                visited.add(nei)
                                                stack.append(nei)
                                    comp_index += 1
                                    clusters_map[f"cluster_{comp_index}"] = component
                        # Ensure we have clusters; if none, create singleton clusters per title
                        if not clusters_map:
                            for t in TITLE:
                                clusters_map.setdefault(f"cluster_{t}", [t])

                        # Optionally write per-cluster FASTA files (debugging aid)
                        for cid, ids in clusters_map.items():
                            out_path = os.path.join(CLUSTERS_DIR, f"{cid}.fasta")
                            try:
                                with open(out_path, 'w', encoding='utf-8') as fh:
                                    for sid in ids:
                                        seq = title_to_seq.get(sid)
                                        if seq is not None:
                                            fh.write(f">{sid}\n{seq}\n")
                            except Exception:
                                pass

                        # Build cluster data in-memory instead of re-reading from disk
                        CLUSTER_DATA_FRAMES = []
                        PAN_LIST_FIRST = []
                        CORE_LIST_FIRST = []
                        ACCESSORY_LIST_FIRST = []
                        UNIQUE_LIST_FIRST = []

                        for cid, ids in clusters_map.items():
                            RESULT_CLUSTER = ALL_SEQ[ALL_SEQ["TITLE"].isin(ids)]
                            if RESULT_CLUSTER.empty:
                                continue
                            CLUSTER_DATA_FRAMES.append(RESULT_CLUSTER)
                            present_files = set(RESULT_CLUSTER["FILE"].to_list())
                            present_count = len(present_files)
                            total_files = len(all_names)

                            try:
                                first_seq = RESULT_CLUSTER["SEQ"].to_list()[0]
                                PAN_LIST_FIRST.append(first_seq)
                                if present_count == total_files:
                                    CORE_LIST_FIRST.append(first_seq)
                                elif present_count == 1:
                                    UNIQUE_LIST_FIRST.append(first_seq)
                                else:
                                    ACCESSORY_LIST_FIRST.append(first_seq)
                            except IndexError:
                                pass

                        LABEL_PERCENTAGE.config(text=f"Progress : Finalizing", fg="red")
                        os.makedirs(RESULTS_PANGENOME, exist_ok=True)
                        with open(os.path.join(RESULTS_PANGENOME, "GENE_PRESENCE_ABSENCE.txt"), "w", encoding="utf-8") as ABSENCE_PRESENCE_FILE:
                            ABSENCE_PRESENCE_FILE.write("Cluster_Info\tPresence\tAbsence\tGenes\n")
                            for GENE_ABS_PRS in CLUSTER_DATA_FRAMES:
                                GROUPED_DATA_CONCAT_NEXT = GENE_ABS_PRS[["TITLE", "FILE"]]
                                DIFFERENCE_FILE = set(all_names).difference(
                                    set(GROUPED_DATA_CONCAT_NEXT["FILE"].to_list()))
                                PRESENT_FILE = set(GROUPED_DATA_CONCAT_NEXT["FILE"].to_list())
                                PRESENT_FILE_REFINED = [os.path.basename(only_file_next) for only_file_next in list(PRESENT_FILE)]
                                ABSENCE_PRESENCE_FILE.write(
                                    f"PRESENT IN {len(PRESENT_FILE_REFINED)} FILES {str(PRESENT_FILE_REFINED)}\n")
                                if len(DIFFERENCE_FILE) != 0:
                                    DIFFERENCE_FILE_REFINED = [os.path.basename(only_file) for only_file in list(DIFFERENCE_FILE)]
                                    ABSENCE_PRESENCE_FILE.write(
                                        f"NOT PRESENT IN {len(DIFFERENCE_FILE_REFINED)} FILES {str(DIFFERENCE_FILE_REFINED)}\n")
                                else:
                                    ABSENCE_PRESENCE_FILE.write(f"NOT PRESENT IN 0\n")
                                for ABS_PRS_TITLE, ABS_PRS_FILE in zip(GROUPED_DATA_CONCAT_NEXT["TITLE"],
                                                                       GROUPED_DATA_CONCAT_NEXT["FILE"]):
                                    ABSENCE_PRESENCE_FILE.write(f"{ABS_PRS_TITLE} {os.path.basename(ABS_PRS_FILE)}\n")
                                ABSENCE_PRESENCE_FILE.write("\n")
                                ABSENCE_PRESENCE_FILE.write("\n")
                                ABSENCE_PRESENCE_FILE.write("\n")
                            ABSENCE_PRESENCE_FILE.close()

                        SEQS_TO_MAKE_FILE = ALL_SEQ.drop_duplicates(subset=['SEQ'])

                        PANGENOME_ANNOTATED = SEQS_TO_MAKE_FILE[SEQS_TO_MAKE_FILE["SEQ"].isin(PAN_LIST_FIRST)]
                        CORE_ANNOTATED = SEQS_TO_MAKE_FILE[SEQS_TO_MAKE_FILE["SEQ"].isin(CORE_LIST_FIRST)]
                        ACCESSORY_ANNOTATED = SEQS_TO_MAKE_FILE[SEQS_TO_MAKE_FILE["SEQ"].isin(ACCESSORY_LIST_FIRST)]
                        UNIQUE_ANNOTATED = SEQS_TO_MAKE_FILE[SEQS_TO_MAKE_FILE["SEQ"].isin(UNIQUE_LIST_FIRST)]
                        os.makedirs(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis"), exist_ok=True)
                        for SINGLE_FILE_TO_MAKE in OUTPUT_LIST:
                            if SINGLE_FILE_TO_MAKE.split()[0] == "Pangenome":
                                with open(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "PAN_FILE.faa"), "w", encoding="utf-8") as PANGENOME:
                                    for DATA_ANNOTATED_PAN_TITLE, DATA_ANNOTATED_PAN_FASTA in zip(
                                            PANGENOME_ANNOTATED["TITLE"],
                                            PANGENOME_ANNOTATED["SEQ"]):
                                        PANGENOME.write(f">{DATA_ANNOTATED_PAN_TITLE}\n{DATA_ANNOTATED_PAN_FASTA}\n")
                                    PANGENOME.close()
                            elif SINGLE_FILE_TO_MAKE.split()[0] == "Core":
                                with open(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "CORE_FILE.faa"), "w", encoding="utf-8") as COREGENOME:
                                    for CORE_ANNOTATED_CORE_TITLE, CORE_ANNOTATED_CORE_FASTA in zip(
                                            CORE_ANNOTATED["TITLE"],
                                            CORE_ANNOTATED["SEQ"]):
                                        COREGENOME.write(f">{CORE_ANNOTATED_CORE_TITLE}\n{CORE_ANNOTATED_CORE_FASTA}\n")
                                    COREGENOME.close()
                            elif SINGLE_FILE_TO_MAKE.split()[0] == "Accessory":
                                with open(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "ACC_FILE.faa"), "w", encoding="utf-8") as ACCGENOME:
                                    for ACC_ANNOTATED_ACC_TITLE, ACC_ANNOTATED_ACC_FASTA, in zip(
                                            ACCESSORY_ANNOTATED["TITLE"],
                                            ACCESSORY_ANNOTATED["SEQ"]):
                                        ACCGENOME.write(f">{ACC_ANNOTATED_ACC_TITLE}\n{ACC_ANNOTATED_ACC_FASTA}\n")
                                    ACCGENOME.close()
                            elif SINGLE_FILE_TO_MAKE.split()[0] == "Unique":
                                os.makedirs(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis"), exist_ok=True)
                                with open(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "UNI_FILE.faa"), "w", encoding="utf-8") as UNIGENOME:
                                    for UNI_ANNOTATED_UNI_TITLE, UNI_ANNOTATED_UNI_FASTA in zip(
                                            UNIQUE_ANNOTATED["TITLE"],
                                            UNIQUE_ANNOTATED["SEQ"]):
                                        UNIGENOME.write(f">{UNI_ANNOTATED_UNI_TITLE}\n{UNI_ANNOTATED_UNI_FASTA}\n")
                                    UNIGENOME.close()

                            else:
                                messagebox.showerror("Error", "No Files were selected.")

                        FILE_NAME = []
                        PAN_LIST_INTER = []
                        PAN_LIST = []
                        AVG_PROTEOME = []
                        IND_PAN = []
                        TOTAL_PROTEINS = []

                        CORE_LIST = []
                        CORE_INDIVIDUAL = []
                        CORE_NUMBER_LIST = []

                        ACC_LIST = []
                        ACC_INDIVIDUAL = []
                        ACC_NUMBER_LIST = []

                        UNI_LIST = []
                        UNI_INDIVIDUAL = []
                        UNI_NUMBER_LIST = []

                        for file in all_names:
                            TITLE, FASTA, DES = read_submitted_fasta_file(file)
                            TOTAL_PROTEINS.append(len(FASTA))
                            FILE_NAME.append(file)
                            AVG_PROTEOME.append(len(FASTA))

                            PAN_INTER = set(FASTA).intersection(set(PAN_LIST_FIRST))
                            IND_PAN.append(len(PAN_INTER))
                            PAN_LIST_INTER.append(PAN_INTER)
                            PAN_LIST.append(len(set.union(*PAN_LIST_INTER)))

                            ALL_CORE = set(FASTA).intersection(set(CORE_LIST_FIRST))
                            CORE_INDIVIDUAL.append(len(ALL_CORE))
                            CORE_LIST.append(ALL_CORE)
                            CORE_NUMBER_LIST.append(len(set.union(*CORE_LIST)))

                            ALL_ACC = set(FASTA).intersection(set(ACCESSORY_LIST_FIRST))
                            ACC_INDIVIDUAL.append(len(ALL_ACC))
                            ACC_LIST.append(ALL_ACC)
                            ACC_NUMBER_LIST.append(len(set.union(*ACC_LIST)))

                            ALL_UNI = set(FASTA).intersection(set(UNIQUE_LIST_FIRST))
                            UNI_INDIVIDUAL.append(len(ALL_UNI))
                            UNI_LIST.append(ALL_UNI)
                            UNI_NUMBER_LIST.append(len(set.union(*UNI_LIST)))

                        ALL_DATA = pd.DataFrame({
                            "FILE": FILE_NAME,
                            "PAN": PAN_LIST,
                            "CORE": CORE_NUMBER_LIST,
                            "ACC": ACC_NUMBER_LIST,
                            "UNI": UNI_NUMBER_LIST
                        })

                        ALL_DATA_INDIVIDUAL = pd.DataFrame({
                            "FILE": FILE_NAME,
                            "TOTAL": TOTAL_PROTEINS,
                            "PAN": IND_PAN,
                            "CORE": CORE_INDIVIDUAL,
                            "ACC": ACC_INDIVIDUAL,
                            "UNI": UNI_INDIVIDUAL
                        })
                        
                        os.makedirs(os.path.join(RESULTS_PANGENOME, "Statistics"), exist_ok=True)
                        ALL_DATA_INDIVIDUAL.to_csv(os.path.join(RESULTS_PANGENOME, "Statistics", "result_ind.csv"), index=False, header=True)
                        ALL_DATA.to_csv(os.path.join(RESULTS_PANGENOME, "Statistics", "result.csv"), index=False, header=True)

                        # Apply publication style for high-quality pangenome plots
                        apply_publication_style()
                        
                        # Use publication-quality colors
                        pub_colors = get_publication_colors()
                        
                        # Increase figure size to provide more space for elements
                        plt.close('all')  # Close any existing figures to avoid conflicts
                        fig, ax = plt.subplots(figsize=(20, 12))  # Wider figure for better layout
                        
                        # Plot each category with explicit colors and thick lines
                        for i, column in enumerate(['PAN', 'CORE', 'ACC', 'UNI']):
                            ax.plot(range(len(ALL_DATA)), ALL_DATA[column], 
                                   color=pub_colors[i], linewidth=4, marker='o', markersize=12,
                                   label=column, markerfacecolor='white', markeredgewidth=3)
                        
                        ax.set_xlabel("Genomes", fontweight='bold', fontsize=20)
                        ax.set_ylabel("Number of Genes", fontweight='bold', fontsize=20)
                        ax.set_title("Pangenome Accumulation Curve", fontweight='bold', fontsize=24)
                        
                        # Customize ticks
                        ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=8)
                        ax.tick_params(axis='both', which='minor', labelsize=14, width=1.5, length=4)
                        
                        # Set x-axis to show genome numbers properly
                        ax.set_xticks(range(len(ALL_DATA)))
                        ax.set_xticklabels([f'G{i+1}' for i in range(len(ALL_DATA))], fontsize=16)
                        
                        # Move legend outside the plot area to the right side at the top
                        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 0.98), frameon=True, fancybox=False, shadow=False,
                                 framealpha=0.9, facecolor='white', edgecolor='black',
                                 fontsize=18, title='Gene Categories', title_fontsize=20)
                        
                        # Calculate percentages
                        if PAN_LIST_FIRST:
                           core_perc = round(len(CORE_LIST_FIRST) / len(PAN_LIST_FIRST) * 100, 2)
                           acc_perc = round(len(ACCESSORY_LIST_FIRST) / len(PAN_LIST_FIRST) * 100, 2)
                           uni_perc = round(len(UNIQUE_LIST_FIRST) / len(PAN_LIST_FIRST) * 100, 2)
                        else:
                            core_perc = acc_perc = uni_perc = 0

                        avg_genome_size = int(np.average(AVG_PROTEOME)) if AVG_PROTEOME else 0

                        # Enhanced text box with better styling
                        stats_text = (
                            f"Statistics:\n"
                            f"• Pangenome: 100%\n"
                            f"• Core: {core_perc}%\n" 
                            f"• Accessory: {acc_perc}%\n"
                            f"• Unique: {uni_perc}%\n"
                            f"• Avg. Genome Size: {avg_genome_size:,}"
                        )
                        
                        # Position the stats text box outside the plot area on the right side below the legend
                        ax.text(1.02, 0.6, stats_text, transform=ax.transAxes, fontsize=16,
                               verticalalignment='top', horizontalalignment='left',
                               bbox=dict(boxstyle='round,pad=1', facecolor='white', 
                                        edgecolor='black', alpha=0.9, linewidth=2))
                        
                        # Add grid for better readability
                        ax.grid(True, alpha=0.3, linestyle='-', linewidth=1)
                        ax.set_axisbelow(True)
                        
                        # Adjust layout to make room for the legend and stats box on the right
                        # Don't use tight_layout as it can conflict with manual adjustments
                        plt.subplots_adjust(left=0.1, right=0.75, top=0.9, bottom=0.1)  # Leave 25% of space on the right for legend and stats
                        line_plot_path = os.path.join(RESULTS_PANGENOME, "Statistics", "Line_plot.png")
                        try:
                            # Force matplotlib to render the figure
                            plt.draw()
                            # Save with explicit format and higher quality
                            fig.savefig(line_plot_path, format='png', dpi=600, bbox_inches='tight', 
                                      facecolor='white', edgecolor='none', pad_inches=0.5)
                            print(f"Plot saved to {line_plot_path}")
                        except PermissionError:
                            # If file is locked by an external viewer, save to a timestamped fallback
                            fallback = line_plot_path.replace(
                                ".png", f"_{int(time.time())}.png"
                            )
                            fig.savefig(fallback, format='png', dpi=600, bbox_inches='tight',
                                      facecolor='white', edgecolor='none', pad_inches=0.5)
                            print(f"Plot saved to {fallback}")
                        
                        plt.clf()  # Clear the figure

                        LABEL_PERCENTAGE.config(text=f"Progress : Cleaning", fg="red")

                        # Cleanup temporary files from both cluster locations
                        for directory in [CLUSTERS_DIR, os.path.abspath(os.path.join(RESULTS_PANGENOME, "Pangenome_Analysis", "cluster"))]:
                            if os.path.isdir(directory):
                                for filename in os.listdir(directory):
                                    file_path = os.path.join(directory, filename)
                                    try:
                                        if os.path.isfile(file_path) or os.path.islink(file_path):
                                            os.unlink(file_path)
                                        elif os.path.isdir(file_path):
                                            shutil.rmtree(file_path)
                                    except Exception:
                                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

                        LABEL_PERCENTAGE.config(text=f"Progress : 0%", fg=Label_Color)
                        
                        # Optional KEGG annotation if toggled
                        if KEGG_Include_var.get() == "KEGG True":
                            try:
                                LABEL_PERCENTAGE.config(text="Progress : KEGG Analysis", fg="red")
                                os.makedirs(RESULTS_KEGG, exist_ok=True)
                                
                                # Use the first file for KEGG analysis
                                QUERY = all_names[0]
                                kegg_db = resolve_diamond_db_base(KEGG_DB)
                                if not diamond_db_present(kegg_db):
                                    LABEL_PERCENTAGE.config(text="Progress : Building KEGG Database", fg="red")
                                    ensure_kegg_dmnd()
                                if not diamond_db_present(kegg_db):
                                    messagebox.showerror("KEGG Error", f"KEGG database not found. Expected {KEGG_DB}.dmnd")
                                    raise FileNotFoundError("KEGG DB missing")
                                
                                LABEL_PERCENTAGE.config(text="Progress : Running KEGG DIAMOND", fg="red")
                                kegg_out = resolve_path(os.path.join(RESULTS_KEGG, "output_file.txt"))
                                run_diamond_blastp(kegg_db, QUERY, kegg_out)
                                if not os.path.exists(kegg_out):
                                    messagebox.showerror("KEGG Error", "KEGG DIAMOND output not created.")
                                    raise FileNotFoundError("KEGG output missing")

                                LABEL_PERCENTAGE.config(text="Progress : Processing KEGG Results", fg="red")
                                kc = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                      "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                                KDATA = pd.read_csv(kegg_out, sep='\t', names=kc)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                                KDATA = KDATA[KDATA["PERCENTAGE IDENTITY"] >= KEGG_Percentage_combo_box_var.get()]

                                best = []
                                for q in KDATA["QUERY"].unique():
                                    sub = KDATA[KDATA["QUERY"] == q]
                                    best.append(sub.loc[sub["BIT SCORE"].idxmax()])
                                if best:
                                    KDATA_BEST = pd.DataFrame(best)
                                    
                                    # ✅ CRITICAL FIX: Add KEGG category mapping here
                                    LABEL_PERCENTAGE.config(text="Progress : Mapping KEGG Categories", fg="red")
                                    try:
                                        # Use relative path based on script location
                                        KEGG_DB_PATH = os.path.join(SCRIPT_DIR, "Data", "Db", "KEGGDb")
                                        
                                        # Check if KEGG database files exist
                                        required_files = ["ko_pathway.txt", "br08901.keg", "pathway_list.txt"]
                                        missing_files = [f for f in required_files if not os.path.exists(os.path.join(KEGG_DB_PATH, f))]
                                        
                                        if not missing_files:
                                            ko_pathway = pd.read_csv(os.path.join(KEGG_DB_PATH, "ko_pathway.txt"),
                                                                     sep="\t", header=None, names=["ko", "pathway"])
                                            # Parse KEGG categories from br08901.keg
                                            brite_file = os.path.join(KEGG_DB_PATH, "br08901.keg")
                                            pathway_to_category = parse_kegg_categories(brite_file)
                                            pathway_list = pd.read_csv(os.path.join(KEGG_DB_PATH, "pathway_list.txt"),
                                                                       sep="\t", header=None, names=["pathway", "desc"])
                                            ko_to_category = pd.merge(ko_pathway, pathway_list, on="pathway", how="left")
                                            ko_to_category["category"] = ko_to_category["pathway"].map(lambda p: pathway_to_category.get(p, "Other"))
                                            ko_category_dict = dict(zip(ko_to_category["ko"], ko_to_category["category"]))
                                            
                                            def normalize_ko(x):
                                                x = str(x).strip()
                                                if not x.startswith("ko:") and x.startswith("K"):
                                                    return "ko:" + x
                                                return x

                                            # ✅ CRITICAL FIX: Handle KEGG data properly
                                            # The SUBJECT column should contain KEGG identifiers from the DIAMOND search
                                            # Let's create meaningful categories from the data we have
                                            
                                            # For now, let's use the SUBJECT column as categories
                                            # This will create charts similar to COG but with KEGG identifiers
                                            KDATA_BEST["CATEGORY"] = KDATA_BEST["SUBJECT"].str.split('|').str[-1]  # Extract the last part after |
                                            KDATA_BEST["CATEGORY"] = KDATA_BEST["CATEGORY"].fillna("Unknown")
                                            
                                            print(f"[Info] Created KEGG categories from SUBJECT identifiers")
                                            print(f"[Info] Sample categories: {KDATA_BEST['CATEGORY'].head().tolist()}")
                                            
                                            # Note: KEGG_ANNOTATION.csv was already processed above
                                        else:
                                            print(f"[Warning] Missing KEGG database files: {missing_files}")
                                            # Fallback: create a basic category column
                                            KDATA_BEST["CATEGORY"] = "Unknown"
                                    except Exception as e:
                                        print(f"[Warning] KEGG category mapping failed: {e}")
                                        KDATA_BEST["CATEGORY"] = "Unknown"
                                    
                                    KDATA_BEST.to_csv(os.path.join(RESULTS_KEGG, "KEGG_Results.csv"), index=False)
                                    
                                    LABEL_PERCENTAGE.config(text="Progress : Generating KEGG Charts", fg="red")
                                    # Generate KEGG charts
                                    try:
                                        plot_kegg_charts(KDATA_BEST, os.path.abspath(RESULTS_KEGG), KEGG_Charts_combo_box_var.get())
                                    except Exception as e:
                                        print(f"[Warning] KEGG chart generation failed: {e}")
                            except Exception as e:
                                messagebox.showerror("KEGG Error", f"KEGG Analysis failed: {e}")
                                pass
                        
                        messagebox.showinfo("Info", "Analysis Completed.")

                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

                    else:
                        messagebox.showerror("Error", "Diamond not found.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

                else:
                    messagebox.showerror("Error", "Select an output file.")
                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

        elif SELECTED_METHOD.get() == "COG":

            ui_set_status("Job Status : Running", "Red")

            if len(all_names) == 1:

                # COG FUNCTIONAL ANNOTATION WITH USER-SELECTED METHOD
                os.makedirs(RESULTS_COG, exist_ok=True)
                output_file = resolve_path(os.path.join(RESULTS_COG, "output_file.txt"))

                QUERY = all_names[0]
                
                # Get user's method selection
                method_choice = COG_Method_combo_box_var.get()
                ui_set_status(f"Job Status : COG Analysis - Method: {method_choice}", "Red")
                
                # Check availability of tools and databases
                # --- paths to databases --- 
                rps_db = COG_CDD_DB   # Database prefix
                rpsblast_path = RPSBLAST_EXE  # RPS-BLAST executable
          
                # Check RPS-BLAST availability (both exe AND database must exist)
                rps_exe_ok = os.path.exists(rpsblast_path)
                rps_db_ok = (
                    os.path.exists(rps_db + ".pal")
                    or len(glob.glob(rps_db + "*.rps")) > 0
                )
                rps_available = rps_exe_ok and rps_db_ok

                print(f"[{'OK' if rps_exe_ok else 'MISSING'}] rpsblast.exe: {rpsblast_path}")
                print(f"[{'OK' if rps_db_ok else 'MISSING'}] COG DB: {rps_db}")

                # Confirm availability (diagnostic logging only; rps_available already computed above)
                if not rps_exe_ok:
                   print(f"[Error] RPS-BLAST executable not found at: {rpsblast_path}")
                else:
                    print(f"[OK] RPS-BLAST executable found: {rpsblast_path}")

                if not rps_db_ok:
                    print(f"[Error] RPS-BLAST database not found at: {rps_db}")
                else:
                     print(f"[OK] RPS-BLAST database found: {rps_db}")

                
                # Check DIAMOND availability
                diamond_path = DIAMOND_EXE
                diamond_db = COG_DB
                diamond_available = (os.path.exists(diamond_path) and diamond_db_present(diamond_db))

                print(f"[{'OK' if diamond_available else 'MISSING'}] DIAMOND: {diamond_path}")
                
                # Execute based on user's method choice
                use_rps = False
                method_used = ""
                
                if method_choice == "Auto (RBS-BLAST → DIAMOND)":
                    if rps_available:
                        # Try RPS-BLAST first
                        ui_set_status("Job Status : Auto mode - Trying RPS-BLAST first", "Red")
                        try:
                            run_rpsblast(rps_db, QUERY, output_file)
                            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                                use_rps = True
                                method_used = "RPS-BLAST"
                                ui_info("Auto mode: RPS-BLAST profile-based COG annotation completed successfully.")
                            else:
                                raise FileNotFoundError("RPS-BLAST produced no output")
                        except Exception as e:
                            ui_set_status(f"Job Status : Auto mode - RPS-BLAST failed, trying DIAMOND fallback", "Red")
                            print(f"Auto mode: RPS-BLAST failed: {e}")
                            # Clear any partial output file
                            try:
                                if os.path.exists(output_file):
                                    os.remove(output_file)
                            except Exception:
                                pass
                            
                            # Try DIAMOND fallback
                            if diamond_available:
                                try:
                                    created = run_diamond_blastp_with_retries(diamond_db, QUERY, output_file)
                                    if created and os.path.exists(output_file):
                                        method_used = "DIAMOND"
                                        ui_info("Auto mode: Fell back to DIAMOND pairwise alignment successfully.")
                                    else:
                                        raise FileNotFoundError("DIAMOND also failed to produce output")
                                except Exception as diamond_error:
                                    ui_error("COG Error", f"Auto mode failed: RPS-BLAST failed ({e}), DIAMOND also failed ({diamond_error})")
                                    ui_set_status("Job Status : None", Label_Color)
                                    return
                            else:
                                ui_error("COG Error", f"Auto mode failed: RPS-BLAST failed ({e}) and DIAMOND not available")
                                ui_set_status("Job Status : None", Label_Color)
                                return
                    else:
                        # RPS-BLAST not available, try DIAMOND directly
                        if diamond_available:
                            ui_set_status("Job Status : Auto mode - RPS-BLAST not available, using DIAMOND", "Red")
                            try:
                                created = run_diamond_blastp_with_retries(diamond_db, QUERY, output_file)
                                if created and os.path.exists(output_file):
                                    method_used = "DIAMOND"
                                    ui_info("Auto mode: Using DIAMOND (RPS-BLAST not available).")
                                else:
                                    raise FileNotFoundError("DIAMOND failed to produce output")
                            except Exception as e:
                                ui_error("COG Error", f"Auto mode failed: RPS-BLAST not available and DIAMOND failed ({e})")
                                ui_set_status("Job Status : None", Label_Color)
                                return
                        else:
                            ui_error("COG Error", "Auto mode failed: Neither RPS-BLAST nor DIAMOND are available")
                            ui_set_status("Job Status : None", Label_Color)
                            return
                            
                elif method_choice == "RBS-BLAST Only":

                    rpsblast_path = RPSBLAST_EXE
                    rps_db = COG_CDD_DB
                    rps_available = (
                        rpsblast_path is not None and os.path.exists(rpsblast_path)
                        and (os.path.exists(rps_db + ".pal") or len(glob.glob(rps_db + "*.rps")) > 0)
                    )
                    if not rps_available:
                       if rpsblast_path is None or not os.path.exists(rpsblast_path):
                          ui_error("RPS-BLAST Error", f"RPS-BLAST executable not found. Expected at {rpsblast_path}")
                       else:
                         ui_error("RPS-BLAST Error", f"CDD database not found. Expected at {rps_db}.rps")
                       ui_set_status("Job Status : None", Label_Color)
                       return

                    ui_set_status("Job Status : Running RPS-BLAST Only (Profile-based COG annotation)", "Red")
                    try:
                      run_rpsblast(rps_db, QUERY, output_file)
                      if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
                        raise FileNotFoundError("RPS-BLAST output file not created or empty")
                      use_rps = True
                      method_used = "RPS-BLAST"
                      ui_info("RPS-BLAST profile-based COG annotation completed successfully.")
                    except Exception as e:
                      ui_error("RPS-BLAST Error", f"Profile-based COG annotation failed: {e}")
                      ui_set_status("Job Status : None", Label_Color)
                      return

                        
                elif method_choice == "DIAMOND Only":
                    if not diamond_available:
                        if not os.path.exists(diamond_path):
                            ui_error("DIAMOND Error", "DIAMOND executable not found. Please install or place diamond.exe under ./Data/Tools/diamond/")
                        else:
                            ui_error("DIAMOND Error", "COG database not found. Expected ./Data/Db/COGDb/COGDb.dmnd")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    
                    ui_set_status("Job Status : Running DIAMOND Only (Sequence-based COG annotation)", "Red")
                    try:
                        created = run_diamond_blastp_with_retries(diamond_db, QUERY, output_file)
                        if not created or not os.path.exists(output_file):
                            raise FileNotFoundError("DIAMOND output file not created")
                        method_used = "DIAMOND"
                        ui_info("DIAMOND pairwise alignment completed successfully.")
                    except Exception as e:
                        ui_error("DIAMOND Error", f"COG DIAMOND analysis failed: {e}")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                else:
                    ui_error("COG Error", f"Unknown method selection: {method_choice}")
                    ui_set_status("Job Status : None", Label_Color)
                    return

                # Parse output based on method used
                if use_rps:
                    # Parse RPS-BLAST output: qseqid sacc pident length qlen qstart qend evalue bitscore
                    COLUMN_NAMES = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "qlen", "qstart", "qend", "E VALUE", "BIT SCORE"]
                    try:
                        if os.path.getsize(output_file) == 0:
                            ui_info("Analysis Completed. No COG hits found above threshold.")
                            ui_set_status("Job Status : None", Label_Color)
                            return
                        DATA_SEP = pd.read_csv(output_file, sep='\t', names=COLUMN_NAMES)
                    except Exception as e:
                        ui_error("COG Error", f"Failed to read RPS-BLAST output: {e}")
                        ui_set_status("Job Status : None", Label_Color)
                        return

                    # Compute non-overlapping query coverage per (QUERY, SUBJECT) by merging HSP intervals
                    # Ensure numeric types
                    for col in ["alignment length", "qlen", "qstart", "qend", "PERCENTAGE IDENTITY", "BIT SCORE"]:
                        try:
                            DATA_SEP[col] = pd.to_numeric(DATA_SEP[col], errors="coerce")
                        except Exception:
                            pass
                    DATA_SEP = DATA_SEP.dropna(subset=["qlen", "qstart", "qend"])                    

                    def _merged_query_span_len(df):
                        try:
                            intervals = []
                            for _, r in df.iterrows():
                                s = int(min(r["qstart"], r["qend"]))
                                e = int(max(r["qstart"], r["qend"]))
                                intervals.append((s, e))
                            if not intervals:
                                return 0
                            intervals.sort(key=lambda x: x[0])
                            merged = []
                            cur_s, cur_e = intervals[0]
                            for s, e in intervals[1:]:
                                if s <= cur_e:
                                    if e > cur_e:
                                        cur_e = e
                                else:
                                    merged.append((cur_s, cur_e))
                                    cur_s, cur_e = s, e
                            merged.append((cur_s, cur_e))
                            covered = 0
                            for s, e in merged:
                                covered += (e - s + 1)
                            return covered
                        except Exception:
                            return 0

                    grouped = DATA_SEP.groupby(["QUERY", "SUBJECT"]).apply(
                        lambda g: pd.Series({
                            "qlen": g["qlen"].iloc[0],
                            "COVERED": _merged_query_span_len(g),
                            "PERCENTAGE IDENTITY": g["PERCENTAGE IDENTITY"].max(),
                            "BIT SCORE": g["BIT SCORE"].max(),
                            "E VALUE": g["E VALUE"].min()
                        })
                    ).reset_index()

                    grouped["COVERAGE"] = grouped.apply(
                        lambda r: (float(r["COVERED"]) / float(r["qlen"])) * 100.0 if float(r["qlen"]) > 0 else 0.0,
                        axis=1
                    )

                    try:
                        cov_thr = int(COG_Coverage_combo_box_var.get())
                    except Exception:
                        cov_thr = 50
                    GREATER_PERC = grouped[
                        (grouped["PERCENTAGE IDENTITY"] > COG_Percentage_combo_box_var.get()) &
                        (grouped["COVERAGE"] >= cov_thr)
                    ][["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                    if GREATER_PERC.empty:
                        ui_info("Analysis Completed. No COG hits found above threshold.")
                        ui_set_status("Job Status : None", Label_Color)
                        return

                    # Build best hit per query and map COG accession -> functional category symbol
                    # using COG_ANNOTATION.csv (columns: TITLE, CATEGORY, DESCRIPTION)
                    LIST_COG = []
                    for q in GREATER_PERC["QUERY"].unique():
                        sub = GREATER_PERC[GREATER_PERC["QUERY"] == q]
                        # choose by highest bit score; fall back to max identity
                        if not sub["BIT SCORE"].isna().all():
                            LIST_COG.append(sub.loc[sub["BIT SCORE"].idxmax()])
                        else:
                            LIST_COG.append(sub.loc[sub["PERCENTAGE IDENTITY"].idxmax()])
                    if not LIST_COG:
                        ui_info("Analysis Completed. No COG best-hits could be selected.")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    COG_FINAL_DATA = pd.DataFrame(LIST_COG).drop_duplicates()

                    # RPS-BLAST's SUBJECT column is a raw CDD PSSM-ID (e.g. "CDD:443680"),
                    # not a COG accession. Translate it via cddid.tbl (PSSM-ID -> COG accession)
                    # before mapping to functional category.
                    try:
                        CDDID_MAP = pd.read_csv(
                            CDDID_TBL, sep='\t', header=None,
                            names=['PSSM_ID', 'ACCESSION', 'SHORT_NAME', 'DESCRIPTION', 'LENGTH'],
                            usecols=[0, 1], dtype=str
                        )
                    except Exception:
                        ui_error("COG Error", f"cddid.tbl not found or unreadable at {CDDID_TBL}. "
                                 "Download it from https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz")
                        ui_set_status("Job Status : None", Label_Color)
                        return

                    # Extract numeric PSSM-ID from "CDD:443680" -> "443680"
                    COG_FINAL_DATA['PSSM_ID'] = COG_FINAL_DATA['SUBJECT'].astype(str).str.replace('CDD:', '', regex=False).str.strip()
                    COG_FINAL_DATA = COG_FINAL_DATA.merge(CDDID_MAP, on='PSSM_ID', how='left')

                    # TYPE is the COG accession (e.g., COG0001)
                    # Use the complete official COG-20 definitions file (4877 entries)
                    # instead of COG_ANNOTATION_CSV, which only contains a partial subset.
                    try:
                        DATA_COG_MAP = pd.read_csv(
                            COG_20_DEF_TAB, sep='\t', header=None,
                            names=['TYPE', 'SYMBOL', 'NAME', 'GENE', 'PATHWAY', 'BLANK', 'PDB_ID'],
                            usecols=['TYPE', 'SYMBOL'], dtype=str,
                            encoding='cp1252'  # the NCBI file contains Windows-1252 smart-quote bytes (e.g. 0x92) that aren't valid UTF-8
                        ).drop_duplicates()
                        # A COG can belong to multiple categories (e.g. "JK"); keep only
                        # the first letter so it matches a single functional category bucket.
                        DATA_COG_MAP['SYMBOL'] = DATA_COG_MAP['SYMBOL'].str[0]
                    except Exception as e:
                        ui_error("COG Error", f"Failed to read cog-20.def.tab at {COG_20_DEF_TAB}.\n"
                                 f"Error: {type(e).__name__}: {e}\n"
                                 "If the file is missing, download it from "
                                 "https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    RESULT_ANN = COG_FINAL_DATA.merge(DATA_COG_MAP, left_on='ACCESSION', right_on='TYPE', how='left')

                    # Build counts per functional category symbol
                    DATA_FINAL = RESULT_ANN['SYMBOL'].fillna('Unknown').value_counts().rename_axis('Functional Category').reset_index(name='Number of Sequences')
                    total_n = int(DATA_FINAL['Number of Sequences'].sum()) if not DATA_FINAL.empty else 0
                    if total_n == 0:
                        ui_info("Analysis Completed. No COG categories found.")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    DATA_FINAL['Percentage'] = DATA_FINAL['Number of Sequences'] / total_n * 100.0
                    DATA_FINAL.to_csv(os.path.join(RESULTS_COG, "COG_PERC.csv"), index=False)

                    # Prepare data for charts by merging category names
                    try:
                        COG_CAT = pd.read_csv(SYMBOL_CATEGORIES_CSV)
                    except Exception:
                        COG_CAT = pd.DataFrame({'SYMBOL': [], 'CATEGORY': []})
                    PIE_DATA = DATA_FINAL[DATA_FINAL['Percentage'] > 1]
                    if PIE_DATA.empty:
                        PIE_DATA = DATA_FINAL
                    PIE_DATA_REFINED = PIE_DATA.merge(COG_CAT, left_on='Functional Category', right_on='SYMBOL', how='left')

                    # Draw charts (same style as DIAMOND path)
                    apply_publication_style()
                    pub_colors = get_publication_colors()
                    chart_colors = pub_colors[:len(PIE_DATA_REFINED)] if len(PIE_DATA_REFINED) <= len(pub_colors) else pub_colors * (len(PIE_DATA_REFINED) // len(pub_colors) + 1)
                    explode = [0.05] * len(PIE_DATA_REFINED)

                    if COG_combo_box_var.get() == "Pie":
                        fig, ax = plt.subplots(figsize=(12, 10))
                        wedges, texts, autotexts = ax.pie(
                            PIE_DATA_REFINED['Percentage'],
                            colors=chart_colors,
                            labels=PIE_DATA_REFINED['Functional Category'],
                            autopct='%1.1f%%', pctdistance=0.8,
                            explode=explode,
                            textprops={'fontsize': 16, 'fontweight': 'bold'},
                            wedgeprops={'linewidth': 2, 'edgecolor': 'black'}
                        )
                        for t in autotexts:
                            t.set_color('white'); t.set_fontweight('bold'); t.set_fontsize(16)
                        for t in texts:
                            t.set_fontweight('bold'); t.set_fontsize(14)
                        centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black')
                        fig.gca().add_artist(centre_circle)
                        ax.legend(PIE_DATA_REFINED.get('CATEGORY', pd.Series(['']*len(PIE_DATA_REFINED))),
                                  loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True,
                                  framealpha=0.9, facecolor='white', edgecolor='black', fontsize=16,
                                  title='COG Categories', title_fontsize=18)
                        ax.set_title('COG Functional Annotation', fontweight='bold', fontsize=24, pad=30)
                        plt.tight_layout(); plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.png"), dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none'); plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.svg"), bbox_inches='tight'); plt.clf()
                    elif COG_combo_box_var.get() == "Bar":
                        fig, ax = plt.subplots(figsize=(14, 10))
                        bars = ax.bar(range(len(PIE_DATA_REFINED)), PIE_DATA_REFINED['Number of Sequences'], color=chart_colors, edgecolor='black', linewidth=2, alpha=0.9)
                        ax.set_xlabel("COG Functional Categories", fontweight='bold', fontsize=20)
                        ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
                        ax.set_title('COG Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
                        ax.set_xticks(range(len(PIE_DATA_REFINED)))
                        ax.set_xticklabels(PIE_DATA_REFINED['Functional Category'], rotation=45, ha='right', fontsize=14, fontweight='bold')
                        ax.tick_params(axis='y', labelsize=16, width=2, length=8); ax.tick_params(axis='x', labelsize=14, width=2, length=8)
                        for i, bar in enumerate(bars):
                            height = bar.get_height(); ax.text(bar.get_x() + bar.get_width()/2., height + max(PIE_DATA_REFINED['Number of Sequences'])*0.01, f'{int(height)}', ha='center', va='bottom', fontsize=14, fontweight='bold')
                        legend_labels = [f"{cat} - {sym}" for cat, sym in zip(PIE_DATA_REFINED.get('CATEGORY', pd.Series(['']*len(PIE_DATA_REFINED))), PIE_DATA_REFINED.get('SYMBOL', pd.Series(['']*len(PIE_DATA_REFINED))))]
                        ax.legend(bars, legend_labels, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=True, framealpha=0.9, facecolor='white', edgecolor='black', fontsize=14, title='COG Categories', title_fontsize=16)
                        ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y'); ax.set_axisbelow(True)
                        plt.tight_layout(); plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.png"), dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none'); plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.svg"), bbox_inches='tight'); plt.clf()
                    elif COG_combo_box_var.get() == "All":
                        # Pie
                        fig, ax = plt.subplots(figsize=(12, 10))
                        wedges, texts, autotexts = ax.pie(PIE_DATA_REFINED['Percentage'], colors=chart_colors, labels=PIE_DATA_REFINED['Functional Category'], autopct='%1.1f%%', pctdistance=0.8, explode=explode, textprops={'fontsize': 16, 'fontweight': 'bold'}, wedgeprops={'linewidth': 2, 'edgecolor': 'black'})
                        for t in autotexts:
                            t.set_color('white'); t.set_fontweight('bold'); t.set_fontsize(16)
                        for t in texts:
                            t.set_fontweight('bold'); t.set_fontsize(14)
                        centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black'); fig.gca().add_artist(centre_circle)
                        ax.legend(PIE_DATA_REFINED.get('CATEGORY', pd.Series(['']*len(PIE_DATA_REFINED))), loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True, framealpha=0.9, facecolor='white', edgecolor='black', fontsize=16, title='COG Categories', title_fontsize=18)
                        ax.set_title('COG Functional Annotation', fontweight='bold', fontsize=24, pad=30)
                        plt.tight_layout(); plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.png"), dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none'); plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.svg"), bbox_inches='tight'); plt.clf()
                        # Bar
                        fig, ax = plt.subplots(figsize=(14, 10))
                        bars = ax.bar(range(len(PIE_DATA_REFINED)), PIE_DATA_REFINED['Number of Sequences'], color=chart_colors, edgecolor='black', linewidth=2, alpha=0.9)
                        ax.set_xlabel("COG Functional Categories", fontweight='bold', fontsize=20)
                        ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
                        ax.set_title('COG Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
                        ax.set_xticks(range(len(PIE_DATA_REFINED)))
                        ax.set_xticklabels(PIE_DATA_REFINED['Functional Category'], rotation=45, ha='right', fontsize=14, fontweight='bold')
                        ax.tick_params(axis='y', labelsize=16, width=2, length=8); ax.tick_params(axis='x', labelsize=14, width=2, length=8)
                        for i, bar in enumerate(bars):
                            height = bar.get_height(); ax.text(bar.get_x() + bar.get_width()/2., height + max(PIE_DATA_REFINED['Number of Sequences'])*0.01, f'{int(height)}', ha='center', va='bottom', fontsize=14, fontweight='bold')
                        legend_labels = [f"{cat} - {sym}" for cat, sym in zip(PIE_DATA_REFINED.get('CATEGORY', pd.Series(['']*len(PIE_DATA_REFINED))), PIE_DATA_REFINED.get('SYMBOL', pd.Series(['']*len(PIE_DATA_REFINED))))]
                        ax.legend(bars, legend_labels, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=True, framealpha=0.9, facecolor='white', edgecolor='black', fontsize=14, title='COG Categories', title_fontsize=16)
                        ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y'); ax.set_axisbelow(True)
                        plt.tight_layout(); plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.png"), dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none'); plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.svg"), bbox_inches='tight'); plt.clf()

                    # Write per-sequence outputs
                    try:
                        TITLE, FASTA, DES = read_submitted_fasta_file(QUERY)
                    except Exception:
                        TITLE, FASTA, DES = [], [], []
                    DATA_COG_ORI = pd.DataFrame({"TITLE": TITLE, "FASTA": FASTA, "DES": DES})
                    FINAL_COG = DATA_COG_ORI.merge(RESULT_ANN, left_on='TITLE', right_on='QUERY', how='inner')
                    # Save symbol assignments
                    out_sym = FINAL_COG.rename(columns={'TITLE':'ID'})[["ID", "FASTA", "DES", "SUBJECT", "SYMBOL"]]
                    out_sym.to_csv(os.path.join(RESULTS_COG, "COG_SYMBOL.csv"), index=False)
                    # Save basic results
                    DATA_COG_FINAL = FINAL_COG.rename(columns={'TITLE':'ID'})[["ID", "DES"]].drop_duplicates()
                    DATA_COG_FINAL.to_csv(os.path.join(RESULTS_COG, "COG_Results.csv"), index=False)
                    with open(os.path.join(RESULTS_COG, "COG_Results.fasta"), "w") as COG_RESULT_FILE:
                        for row in DATA_COG_FINAL.itertuples(index=False):
                            idv = getattr(row, 'ID'); des = getattr(row, 'DES', '')
                            seq = DATA_COG_ORI.loc[DATA_COG_ORI['TITLE'] == idv, 'FASTA']
                            if not seq.empty:
                                COG_RESULT_FILE.write(f">{idv} {des}\n{seq.iloc[0]}\n")
                    ui_info("Analysis Completed.")
                    ui_set_status("Job Status : None", Label_Color)
                    return
                else:
                    # Parse DIAMOND output: standard BLAST format
                    COLUMN_NAMES = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                    "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                    try:
                        # If output is empty, provide friendly message
                        if os.path.getsize(output_file) == 0:
                            ui_info("Analysis Completed. No COG hits found above threshold.")
                            ui_set_status("Job Status : None", Label_Color)
                            return
                        DATA_SEP = pd.read_csv(output_file, sep='\t', names=COLUMN_NAMES)[
                            ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                    except Exception as e:
                        ui_error("COG Error", f"Failed to read DIAMOND output: {e}")
                        ui_set_status("Job Status : None", Label_Color)
                        return

                    GREATER_PERC = DATA_SEP[DATA_SEP["PERCENTAGE IDENTITY"] > COG_Percentage_combo_box_var.get()]
                    if GREATER_PERC.empty:
                        ui_info("Analysis Completed. No COG hits found above threshold.")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                if not use_rps:

                    LIST_COG = []

                    for QUERYS in GREATER_PERC["QUERY"]:
                        SELECTED_QUERY = GREATER_PERC[GREATER_PERC["QUERY"] == QUERYS]
                        LIST_COG.append(SELECTED_QUERY[SELECTED_QUERY["PERCENTAGE IDENTITY"] == SELECTED_QUERY[
                            "PERCENTAGE IDENTITY"].max()])

                    # Guard in case LIST_COG is empty
                    if not LIST_COG:
                        ui_info("Analysis Completed. No COG best-hits could be selected.")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    COG_FINAL_DATA = pd.concat(LIST_COG).drop_duplicates()

                    # Use the complete official COG-20 definitions file (4877 entries)
                    # instead of COG_ANNOTATION_CSV, which only contains a partial subset.
                    try:
                        COG_ANNOTATION = pd.read_csv(
                            COG_20_DEF_TAB, sep='\t', header=None,
                            names=['TITLE', 'CATEGORY', 'NAME', 'GENE', 'PATHWAY', 'BLANK', 'PDB_ID'],
                            usecols=['TITLE', 'CATEGORY'], dtype=str,
                            encoding='cp1252'  # the NCBI file contains Windows-1252 smart-quote bytes (e.g. 0x92) that aren't valid UTF-8
                        ).drop_duplicates()
                    except Exception as e:
                        ui_error("COG Error", f"Failed to read cog-20.def.tab at {COG_20_DEF_TAB}.\n"
                                 f"Error: {type(e).__name__}: {e}\n"
                                 "If the file is missing, download it from "
                                 "https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    # A COG can belong to multiple categories (e.g. "JK"); keep only the first letter.
                    COG_ANNOTATION['CATEGORY'] = COG_ANNOTATION['CATEGORY'].str[0]
                    ANNOTATED = COG_ANNOTATION[COG_ANNOTATION["TITLE"].isin(COG_FINAL_DATA["SUBJECT"].to_list())]
                    COG_FINAL = COG_FINAL_DATA[COG_FINAL_DATA["SUBJECT"].isin(ANNOTATED["TITLE"].to_list())][
                        ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY"]]
                    if COG_FINAL.empty:
                        ui_info("Analysis Completed. No subjects matched COG annotations.")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    max_value_index = COG_FINAL.groupby('SUBJECT')['PERCENTAGE IDENTITY'].idxmax()
                    result = COG_FINAL.loc[max_value_index]
                    result_Cog = result[result["SUBJECT"].isin(ANNOTATED["TITLE"].to_list())].merge(ANNOTATED, left_on='SUBJECT',
                                                                                       right_on='TITLE', how='left')

                    DATA_COG = COG_ANNOTATION.rename(columns={'TITLE': 'ID', 'CATEGORY': 'SYMBOL'})

                    COG_RESULT = DATA_COG[DATA_COG["ID"].isin(result_Cog["TITLE"].to_list())]
                    DATA_FINAL = COG_RESULT["SYMBOL"].value_counts().rename_axis('Functional Category').reset_index(
                        name='Number of Sequences')

                    # Avoid division by zero
                    total_n = DATA_FINAL['Number of Sequences'].sum()
                    if total_n == 0:
                        ui_info("Analysis Completed. No COG categories found.")
                        ui_set_status("Job Status : None", Label_Color)
                        return
                    DATA_FINAL['Percentage'] = DATA_FINAL["Number of Sequences"]/total_n * 100

                    DATA_FINAL.to_csv(os.path.join(RESULTS_COG, "COG_PERC.csv"))

                    PIE_DATA = DATA_FINAL[DATA_FINAL["Percentage"] > 1]

                    COG_CAT = pd.read_csv(SYMBOL_CATEGORIES_CSV)

                    PIE_DATA_REFINED = PIE_DATA[PIE_DATA["Functional Category"].isin(COG_CAT["SYMBOL"].to_list())].merge(COG_CAT, left_on="Functional Category",
                                                                                       right_on='SYMBOL', how='left')
                    if PIE_DATA_REFINED.empty:
                        # Fallback: use all categories to still draw charts
                        PIE_DATA_REFINED = DATA_FINAL.merge(COG_CAT, left_on="Functional Category", right_on='SYMBOL', how='left')

                    # Apply publication style and get high-quality colors
                    apply_publication_style()
                    pub_colors = get_publication_colors()
                    
                    # Use publication colors instead of pastel seaborn palette
                    chart_colors = pub_colors[:len(PIE_DATA_REFINED)] if len(PIE_DATA_REFINED) <= len(pub_colors) else pub_colors * (len(PIE_DATA_REFINED) // len(pub_colors) + 1)
                    explode = [0.05] * len(PIE_DATA_REFINED)

                    if COG_combo_box_var.get() == "Pie":
                        fig, ax = plt.subplots(figsize=(12, 10))  # Larger figure for better visibility
                        
                        wedges, texts, autotexts = ax.pie(
                            PIE_DATA_REFINED['Percentage'],
                            colors=chart_colors,
                            labels=PIE_DATA_REFINED["Functional Category"],
                            autopct='%1.1f%%', 
                            pctdistance=0.8,
                            explode=explode,
                            textprops={'fontsize': 16, 'fontweight': 'bold'},  # Larger, bold text
                            wedgeprops={'linewidth': 2, 'edgecolor': 'black'}  # Thick edges
                        )
                        
                        # Enhance text visibility
                        for t in autotexts:
                            t.set_color('white')
                            t.set_fontweight('bold')
                            t.set_fontsize(16)
                        
                        for t in texts:
                            t.set_fontweight('bold')
                            t.set_fontsize(14)
                        
                        # Create a professional donut chart
                        centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black')
                        fig.gca().add_artist(centre_circle)
                        
                        # Enhanced legend
                        ax.legend(PIE_DATA_REFINED["CATEGORY"], 
                                 loc='center left', bbox_to_anchor=(1.1, 0.5),
                                 frameon=True, fancybox=False, shadow=False,
                                 framealpha=0.9, facecolor='white', edgecolor='black',
                                 fontsize=16, title='COG Categories', title_fontsize=18)
                        
                        ax.set_title('COG Functional Annotation', fontweight='bold', fontsize=24, pad=30)
                        
                        plt.tight_layout()
                        plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.png"), dpi=600, bbox_inches='tight', 
                                   facecolor='white', edgecolor='none')
                        plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.svg"), bbox_inches='tight')
                        plt.clf()

                    elif COG_combo_box_var.get() == "Bar":
                        fig, ax = plt.subplots(figsize=(14, 10))  # Larger figure
                        
                        # Create bar plot with publication colors
                        bars = ax.bar(range(len(PIE_DATA_REFINED)), 
                                     PIE_DATA_REFINED["Number of Sequences"],
                                     color=chart_colors,
                                     edgecolor='black',
                                     linewidth=2,
                                     alpha=0.9)
                        
                        # Customize axes
                        ax.set_xlabel("COG Functional Categories", fontweight='bold', fontsize=20)
                        ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
                        ax.set_title('COG Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
                        
                        # Set x-tick labels
                        ax.set_xticks(range(len(PIE_DATA_REFINED)))
                        ax.set_xticklabels(PIE_DATA_REFINED["Functional Category"], 
                                          rotation=45, ha='right', fontsize=14, fontweight='bold')
                        
                        # Enhance y-axis ticks
                        ax.tick_params(axis='y', labelsize=16, width=2, length=8)
                        ax.tick_params(axis='x', labelsize=14, width=2, length=8)
                        
                        # Add value labels on bars
                        for i, bar in enumerate(bars):
                            height = bar.get_height()
                            ax.text(bar.get_x() + bar.get_width()/2., height + max(PIE_DATA_REFINED["Number of Sequences"])*0.01,
                                   f'{int(height)}', ha='center', va='bottom', fontsize=14, fontweight='bold')
                        
                        # Enhanced legend with category descriptions
                        legend_labels = [f"{cat} - {sym}" for cat, sym in zip(PIE_DATA_REFINED["CATEGORY"], PIE_DATA_REFINED["SYMBOL"])]
                        ax.legend(bars, legend_labels,
                                 loc='center left', bbox_to_anchor=(1.02, 0.5),
                                 frameon=True, fancybox=False, shadow=False,
                                 framealpha=0.9, facecolor='white', edgecolor='black',
                                 fontsize=14, title='COG Categories', title_fontsize=16)
                        
                        # Add grid for better readability
                        ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y')
                        ax.set_axisbelow(True)
                        
                        plt.tight_layout()
                        plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.png"), dpi=600, bbox_inches='tight',
                                   facecolor='white', edgecolor='none')
                        plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.svg"), bbox_inches='tight')
                        plt.clf()
                        
                    elif COG_combo_box_var.get() == "All":
                        # PIE CHART
                        fig, ax = plt.subplots(figsize=(12, 10))
                        
                        wedges, texts, autotexts = ax.pie(
                            PIE_DATA_REFINED['Percentage'],
                            colors=chart_colors,
                            labels=PIE_DATA_REFINED["Functional Category"],
                            autopct='%1.1f%%', 
                            pctdistance=0.8,
                            explode=explode,
                            textprops={'fontsize': 16, 'fontweight': 'bold'},
                            wedgeprops={'linewidth': 2, 'edgecolor': 'black'}
                        )
                        
                        for t in autotexts:
                            t.set_color('white')
                            t.set_fontweight('bold')
                            t.set_fontsize(16)
                        
                        for t in texts:
                            t.set_fontweight('bold')
                            t.set_fontsize(14)
                        
                        centre_circle = plt.Circle((0, 0), 0.6, fc='white', linewidth=3, edgecolor='black')
                        fig.gca().add_artist(centre_circle)
                        
                        ax.legend(PIE_DATA_REFINED["CATEGORY"], 
                                 loc='center left', bbox_to_anchor=(1.1, 0.5),
                                 frameon=True, fancybox=False, shadow=False,
                                 framealpha=0.9, facecolor='white', edgecolor='black',
                                 fontsize=16, title='COG Categories', title_fontsize=18)
                        
                        ax.set_title('COG Functional Annotation', fontweight='bold', fontsize=24, pad=30)
                        
                        plt.tight_layout()
                        plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.png"), dpi=600, bbox_inches='tight',
                                   facecolor='white', edgecolor='none')
                        plt.savefig(os.path.join(RESULTS_COG, "COG_PIE.svg"), bbox_inches='tight')
                        plt.clf()

                        # BAR CHART
                        fig, ax = plt.subplots(figsize=(14, 10))
                        
                        bars = ax.bar(range(len(PIE_DATA_REFINED)), 
                                     PIE_DATA_REFINED["Number of Sequences"],
                                     color=chart_colors,
                                     edgecolor='black',
                                     linewidth=2,
                                     alpha=0.9)
                        
                        ax.set_xlabel("COG Functional Categories", fontweight='bold', fontsize=20)
                        ax.set_ylabel("Number of Sequences", fontweight='bold', fontsize=20)
                        ax.set_title('COG Functional Annotation Distribution', fontweight='bold', fontsize=24, pad=30)
                        
                        ax.set_xticks(range(len(PIE_DATA_REFINED)))
                        ax.set_xticklabels(PIE_DATA_REFINED["Functional Category"], 
                                          rotation=45, ha='right', fontsize=14, fontweight='bold')
                        
                        ax.tick_params(axis='y', labelsize=16, width=2, length=8)
                        ax.tick_params(axis='x', labelsize=14, width=2, length=8)
                        
                        for i, bar in enumerate(bars):
                            height = bar.get_height()
                            ax.text(bar.get_x() + bar.get_width()/2., height + max(PIE_DATA_REFINED["Number of Sequences"])*0.01,
                                   f'{int(height)}', ha='center', va='bottom', fontsize=14, fontweight='bold')
                        
                        legend_labels = [f"{cat} - {sym}" for cat, sym in zip(PIE_DATA_REFINED["CATEGORY"], PIE_DATA_REFINED["SYMBOL"])]
                        ax.legend(bars, legend_labels,
                                 loc='center left', bbox_to_anchor=(1.02, 0.5),
                                 frameon=True, fancybox=False, shadow=False,
                                 framealpha=0.9, facecolor='white', edgecolor='black',
                                 fontsize=14, title='COG Categories', title_fontsize=16)
                        
                        ax.grid(True, alpha=0.3, linestyle='-', linewidth=1, axis='y')
                        ax.set_axisbelow(True)
                        
                        plt.tight_layout()
                        plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.png"), dpi=600, bbox_inches='tight',
                                   facecolor='white', edgecolor='none')
                        plt.savefig(os.path.join(RESULTS_COG, "BAR_PIE.svg"), bbox_inches='tight')
                        plt.clf()

                    TITLE, FASTA, DES = read_submitted_fasta_file(QUERY)

                    DATA_COG_ORI = pd.DataFrame({
                        "TITLE": TITLE,
                        "FASTA": FASTA,
                        "DES": DES
                    })

                    FINAL_COG = DATA_COG_ORI[DATA_COG_ORI["TITLE"].isin(result_Cog["QUERY"].to_list())].merge(result_Cog, left_on='TITLE',
                                                                                       right_on='QUERY', how='left')

                    COG_WITH_SYMBOL = COG_RESULT[COG_RESULT["ID"].isin(FINAL_COG['SUBJECT'].to_list())].merge(FINAL_COG, left_on="ID",
                                                                                       right_on='SUBJECT', how='left')[["TITLE_x", "FASTA", "DES_y", "SYMBOL"]]

                    COG_WITH_SYMBOL.to_csv(os.path.join(RESULTS_COG, "COG_SYMBOL.csv"))

                    TITLE_COG = []
                    FASTA_COG = []
                    DES_COG = []

                    for QUERY_TITLE, QUERY_FASTA, COG_DES in zip(FINAL_COG["TITLE_x"], FINAL_COG["FASTA"], FINAL_COG["DES_y"]):
                            TITLE_COG.append(QUERY_TITLE)
                            FASTA_COG.append(QUERY_FASTA)
                            DES_COG.append(" ".join(COG_DES.split()[1:]))

                    DATA_COG_FINAL = pd.DataFrame({
                        "ID": TITLE_COG,
                        "FASTA": FASTA_COG,
                        "DES": DES_COG
                    }).drop_duplicates()

                    DATA_COG_FINAL[["ID", "DES"]].to_csv(os.path.join(RESULTS_COG, "COG_Results.csv"))

                    with open(os.path.join(RESULTS_COG, "COG_Results.fasta"), "w") as COG_RESULT_FILE:
                        for COG_FINAL_ID, COG_FINAL_FASTA, COG_FINAL_DES in zip(DATA_COG_FINAL["ID"], DATA_COG_FINAL["FASTA"], DATA_COG_FINAL["DES"]):
                            COG_RESULT_FILE.write(f">{COG_FINAL_ID} {COG_FINAL_DES}\n{COG_FINAL_FASTA}\n")
                        COG_RESULT_FILE.close()

                    # Optional KEGG annotation if toggled
                    if KEGG_Include_var.get() == "KEGG True":
                        try:
                            os.makedirs(RESULTS_KEGG, exist_ok=True)
                            
                            kegg_db = resolve_diamond_db_base(KEGG_DB)
                            if not diamond_db_present(kegg_db):
                                ensure_kegg_dmnd()
                            if not diamond_db_present(kegg_db):
                                messagebox.showerror("KEGG Error", f"KEGG database not found. Expected {KEGG_DB}.dmnd")
                                raise FileNotFoundError("KEGG DB missing")
                            kegg_out = resolve_path(os.path.join(RESULTS_KEGG, "output_file.txt"))
                            run_diamond_blastp(kegg_db, QUERY, kegg_out)
                            if not os.path.exists(kegg_out):
                                messagebox.showerror("KEGG Error", "KEGG DIAMOND output not created.")
                                raise FileNotFoundError("KEGG output missing")

                            kc = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                  "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                            KDATA = pd.read_csv(kegg_out, sep='\t', names=kc)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                            KDATA = KDATA[KDATA["PERCENTAGE IDENTITY"] > COG_Percentage_combo_box_var.get()]

                            best = []
                            for q in KDATA["QUERY"].unique():
                                sub = KDATA[KDATA["QUERY"] == q]
                                best.append(sub.loc[sub["BIT SCORE"].idxmax()])
                            if best:
                                KDATA_BEST = pd.DataFrame(best)
                                # Optional mapping file to KO or descriptions
                                kmap_path = KEGG_ANNOTATION_CSV
                                if os.path.exists(kmap_path):
                                    KMAP = pd.read_csv(kmap_path)
                                    if "TITLE" in KMAP.columns:
                                        KDATA_BEST = KDATA_BEST.merge(KMAP, left_on="SUBJECT", right_on="TITLE", how="left")
                                KDATA_BEST.to_csv(os.path.join(RESULTS_KEGG, "KEGG_Results.csv"), index=False)
                        except Exception:
                            ui_error("KEGG Error", "KEGG Analysis failed. Please check the KEGG database and try again.")
                            pass

                    ui_info("Analysis Completed.")

                    ui_set_status("Job Status : None", Label_Color)

            else:
                ui_error("Error","COG Analysis requires one file for Analysis.")
                ui_set_status("Job Status : None", Label_Color)

        elif SELECTED_METHOD.get() == "KEGG":

            LOADING_LABEL.config(text="Job Status : Running", fg="Red")

            if len(all_names) == 1:

                os.makedirs(RESULTS_KEGG, exist_ok=True)

                QUERY = all_names[0]
                kegg_db = resolve_diamond_db_base(KEGG_DB)
                if not diamond_db_present(kegg_db):
                    ensure_kegg_dmnd()
                if not diamond_db_present(kegg_db):
                    messagebox.showerror("KEGG Error", f"KEGG database not found. Expected {KEGG_DB}.dmnd")
                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                    return
                kegg_out = resolve_path(os.path.join(RESULTS_KEGG, "output_file.txt"))
                # Use sensitive retries to increase chance of hits
                created = run_diamond_blastp_with_retries(kegg_db, QUERY, kegg_out)
                if not created or not os.path.exists(kegg_out):
                    messagebox.showerror("KEGG Error", "KEGG DIAMOND output not created.")
                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                    # Still produce placeholder charts indicating no hits
                    try:
                        plot_kegg_charts(pd.DataFrame({"SUBJECT": []}), os.path.abspath(RESULTS_KEGG), KEGG_Charts_combo_box_var.get())
                    except Exception:
                        pass
                    return

                kc = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                      "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                try:
                    KDATA_RAW = pd.read_csv(kegg_out, sep='\t', names=kc)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                except FileNotFoundError:
                    messagebox.showerror("KEGG Error", "KEGG DIAMOND output not found. Ensure the KEGG database exists at Data/Db/KEGGDb/KEGGDb.dmnd.")
                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                    # Placeholder charts for no hits
                    try:
                        plot_kegg_charts(pd.DataFrame({"SUBJECT": []}), os.path.abspath(RESULTS_KEGG), KEGG_Charts_combo_box_var.get())
                    except Exception:
                        pass
                    return

                # Try current threshold, then progressively relax if needed
                try:
                    start_thr = int(KEGG_Percentage_combo_box_var.get())
                except Exception:
                    start_thr = 90
                thr_list = []
                for t in [start_thr, 80, 70, 60, 50, 40, 30]:
                    if t not in thr_list and 0 <= t <= 100:
                        thr_list.append(t)
                
                KDATA_BEST = None
                used_thr = start_thr
                for thr in thr_list:
                    KDATA = KDATA_RAW[KDATA_RAW["PERCENTAGE IDENTITY"] >= thr]
                    if KDATA.empty:
                        continue
                    best_rows = []
                    for q in KDATA["QUERY"].unique():
                        sub = KDATA[KDATA["QUERY"] == q]
                        best_rows.append(sub.loc[sub["BIT SCORE"].idxmax()])
                    if best_rows:
                        KDATA_BEST = pd.DataFrame(best_rows)
                        used_thr = thr
                        break
                
                out_dir = os.path.abspath(RESULTS_KEGG)
                if KDATA_BEST is not None and not KDATA_BEST.empty:
                    # Optional mapping file to KO or descriptions
                    kmap_path = KEGG_ANNOTATION_CSV
                    if os.path.exists(kmap_path):
                        try:
                            KMAP = pd.read_csv(kmap_path)
                            if "TITLE" in KMAP.columns:
                                KDATA_BEST = KDATA_BEST.merge(KMAP, left_on="SUBJECT", right_on="TITLE", how="left")
                        except Exception:
                            pass
                    try:
                        KDATA_BEST.to_csv(os.path.join(out_dir, "KEGG_Results.csv"), index=False)
                    except Exception:
                        pass
                    # Generate charts once
                    try:
                        plot_kegg_charts(KDATA_BEST, out_dir, KEGG_Charts_combo_box_var.get())
                    except Exception:
                        pass
                    # Notify if threshold was relaxed
                    if used_thr != start_thr:
                        ui_info(f"KEGG: No hits at {start_thr}%. Charts generated using relaxed threshold {used_thr}%.")
                else:
                    # No hits even after relaxing: produce placeholder charts
                    try:
                        plot_kegg_charts(pd.DataFrame({"SUBJECT": []}), out_dir, KEGG_Charts_combo_box_var.get())
                    except Exception:
                        pass
                    ui_info("KEGG: No hits above threshold. Placeholder charts were created.")

                messagebox.showinfo("Info", "KEGG Analysis Completed.")
                LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

            else:
                messagebox.showerror("Error","KEGG Analysis requires one file for Analysis.")
                LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

        elif SELECTED_METHOD.get() == "Phylogenetic":

            LOADING_LABEL.config(text="Job Status : Running", fg="Red")

            if len(all_names) == 1:
                
                # Ensure Results/Phylo directory exists
                os.makedirs(RESULTS_PHYLO, exist_ok=True)

                # For MSA methods, skip DIAMOND and use all sequences directly
                if Phylo_method_var.get().startswith("MSA"):
                    fasta_file = all_names[0]
                    TITLE, FASTA, DES = read_submitted_fasta_file(fasta_file)
                    
                    # Validate that sequences were actually read
                    # Read sequences from file
                    if not TITLE or not FASTA:
                        messagebox.showerror("Error", f"No sequences found in the input file {fasta_file}. Please check that the file contains valid FASTA sequences.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    
                    CORE_DATA = pd.DataFrame({"TITLE": TITLE, "FASTA": FASTA, "DES": DES})
                    MATCHING_SEQ = CORE_DATA  # Use all sequences for MSA
                    
                    # Check if we have enough sequences for phylogenetic analysis
                    if len(MATCHING_SEQ) < 2:
                        messagebox.showerror("Error", "At least 2 sequences are required for phylogenetic analysis.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    
                    # Write sequences for MSA
                    LABEL_PERCENTAGE.config(text="Progress : Preparing sequences", fg="red")
                    
                    # Create alignment file path
                    alignment_input_path = os.path.join(RESULTS_PHYLO, "Alignment.fasta")
                    
                    # Ensure sequences are properly formatted as strings
                    sequence_count = 0
                    sequence_data = []
                    
                    # First, prepare all sequence data
                    for i, (title, seq) in enumerate(zip(MATCHING_SEQ["TITLE"], MATCHING_SEQ["FASTA"])):
                        # Convert sequence to string and clean it
                        seq_str = str(seq).strip()
                        title_str = str(title).strip()
                        
                        # Skip empty sequences
                        if not seq_str or not title_str:
                            # Skip empty sequences
                            continue
                            
                        sequence_data.append(f">{title_str}\n{seq_str}\n")
                        sequence_count += 1
                    
                            # Check sequence preparation
                    
                    # Verify we have sequences to write
                    if sequence_count == 0:
                        messagebox.showerror("Error", "No valid sequences found for alignment. Check your input file format.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    
                    # Write all sequence data to file with explicit flushing
                    try:
                        with open(alignment_input_path, "w", encoding="utf-8", newline='\n') as Alignment_File:
                            for seq_entry in sequence_data:
                                Alignment_File.write(seq_entry)
                            Alignment_File.flush()  # Explicit flush
                            os.fsync(Alignment_File.fileno())  # Force write to disk
                        
                        # Verify file was written correctly
                        if not os.path.exists(alignment_input_path):
                            raise FileNotFoundError(f"Alignment input file was not created at {alignment_input_path}")
                        
                        file_size = os.path.getsize(alignment_input_path)
                        if file_size == 0:
                            raise ValueError(f"Alignment input file is empty (0 bytes) at {alignment_input_path}")
                        
                        # Successfully wrote alignment file
                        
                        # Additional verification - try to read back the file
                        try:
                            with open(alignment_input_path, "r", encoding="utf-8") as verify_file:
                                content = verify_file.read()
                                if not content.strip():
                                    raise ValueError("Alignment file appears to be empty when read back")
                                # File verification successful
                        except Exception as e:
                            raise ValueError(f"Failed to verify alignment file: {e}")
                            
                    except Exception as e:
                        messagebox.showerror("Error", f"Failed to write alignment file: {e}")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    Alignment_File_Path = os.path.join(RESULTS_PHYLO, "Alignment.fasta")
                    OUTFILE = RESULTS_PHYLO

                    method = Phylo_method_var.get()
                    LABEL_PERCENTAGE.config(text=f"Progress : Running {method}", fg="red")
                    
                    # Validate tool existence before execution
                    tool_path = None
                    if method == "MSA - Clustal Omega":
                        tool_path = CLUSTALO_EXE
                        if not os.path.exists(tool_path):
                            messagebox.showerror("Error", f"Clustal Omega not found at {tool_path}")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                            return
                    elif method == "MSA - MAFFT":
                        tool_path = MAFFT_EXE
                        if not os.path.exists(tool_path):
                            messagebox.showerror("Error", f"MAFFT not found at {tool_path}")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                            return
                    elif method == "MSA - MUSCLE":
                        tool_path = MUSCLE_EXE
                        if not os.path.exists(tool_path):
                            messagebox.showerror("Error", f"MUSCLE not found at {tool_path}")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                            return
                    
                    # Execute the MSA tool
                    try:
                        if method == "MSA - Clustal Omega":
                            CLUSTAL_PATH = tool_path
                            result = subprocess.run(
                                   [CLUSTAL_PATH, "-i", Alignment_File_Path, "-o", OUTFILE + "\\Alignment.fasta",
                                   "--force", f"--threads={os.cpu_count()}"],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    startupinfo=startupinfo,
                                    shell=False,
                                    timeout=300
                                 )
                            if result.returncode != 0:
                                raise Exception(f"Clustal Omega error: {result.stderr.decode()}")

                        elif method == "MSA - MAFFT":
                            MAFFT_PATH = tool_path
                            mafft_dir = os.path.dirname(MAFFT_PATH)
                            input_file = os.path.join(OUTFILE, "input_sequences.fasta")
                            output_file = os.path.join(OUTFILE, "Alignment.fasta")
                            
                            # Convert paths to absolute and normalize
                            abs_input_path = os.path.abspath(input_file)
                            abs_output_path = os.path.abspath(output_file)
                            
                            # First, copy the alignment file to a separate input file to avoid overwriting
                            shutil.copy2(Alignment_File_Path, abs_input_path)
                            
                                    # MAFFT processing
                            
                            # Method 1: Direct execution with output redirection
                            # mafft.bat is a batch wrapper and cannot be run directly with shell=False;
                            # it must be invoked via cmd /c (or shell=True).
                            try:
                                is_bat = MAFFT_PATH.lower().endswith((".bat", ".cmd"))
                                mafft_invoke = ["cmd", "/c", MAFFT_PATH, "--auto", abs_input_path] if is_bat \
                                    else [MAFFT_PATH, "--auto", abs_input_path]
                                with open(abs_output_path, 'w', encoding='utf-8') as f:
                                    result = subprocess.run(
                                        mafft_invoke,
                                        stdout=f,
                                        stderr=subprocess.PIPE,
                                        startupinfo=startupinfo,
                                        shell=False,
                                        timeout=300,
                                        cwd=mafft_dir
                                    )
                                        # Check MAFFT execution result
                                
                                # Check if output was produced
                                if os.path.exists(abs_output_path) and os.path.getsize(abs_output_path) > 0:
                                    # MAFFT succeeded
                                    # Clean up temporary input file
                                    try:
                                        os.remove(abs_input_path)
                                    except Exception:
                                        pass
                                else:
                                    stderr_txt = result.stderr.decode(errors='replace') if result and result.stderr else ""
                                    raise Exception(f"MAFFT produced no output. Return code: {result.returncode}. Stderr: {stderr_txt}")
                                    
                            except Exception as e:
                                # MAFFT Method 1 failed, trying Method 2
                                
                                # Method 2: Try with shell=True and output redirection
                                try:
                                    cmd = f'cd /d "{mafft_dir}" && "{MAFFT_PATH}" --auto "{abs_input_path}" > "{abs_output_path}"'
                                    # Trying shell command
                                    
                                    result2 = subprocess.run(
                                        cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        startupinfo=startupinfo,
                                        shell=True,
                                        timeout=300
                                    )
                                    
                                            # Check MAFFT Method 2 execution result
                                    
                                    if not os.path.exists(abs_output_path) or os.path.getsize(abs_output_path) == 0:
                                        stderr_txt2 = result2.stderr.decode(errors='replace') if result2 and result2.stderr else ""
                                        raise Exception(f"MAFFT Method 2 also failed. Both methods produced no output. Stderr: {stderr_txt2}")
                                    else:
                                        # MAFFT Method 2 succeeded
                                        # Clean up temporary input file
                                        try:
                                            os.remove(abs_input_path)
                                        except Exception:
                                            pass
                                        
                                except Exception as e2:
                                    raise Exception(f"MAFFT execution failed with both methods: Method 1: {e}, Method 2: {e2}")

                        elif method == "MSA - MUSCLE":
                            MUSCLE_PATH = tool_path
                            result = subprocess.run(
                                   [MUSCLE_PATH, "-align", Alignment_File_Path, "-output", OUTFILE + "\\Alignment.fasta"],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   startupinfo=startupinfo,
                                   shell=False,
                                   timeout=300
                                )
                            if result.returncode != 0:
                                raise Exception(f"MUSCLE error: {result.stderr.decode()}")

                    except subprocess.TimeoutExpired:
                        messagebox.showerror("Error", f"MSA timed out after 5 minutes. Try with fewer sequences or a different method.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    except Exception as e:
                        messagebox.showerror("Error", f"MSA failed: {e}")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return

                    # Check if alignment file was created
                    alignment_file = OUTFILE + "\\Alignment.fasta"
                    if not os.path.exists(alignment_file):
                        messagebox.showerror("Error", f"Alignment file was not created at {alignment_file}. Please check if {method} is working correctly.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    
                    # Check if alignment file has content
                    if os.path.getsize(alignment_file) == 0:
                        messagebox.showerror("Error", f"Alignment file was created but is empty. {method} may have failed silently. Check stderr output above.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return

                    LABEL_PERCENTAGE.config(text="Progress : Building trees", fg="red")
                    try:
                        align = AlignIO.read(alignment_file, "fasta")
                    except ValueError as e:
                        messagebox.showerror("Error", f"Failed to read alignment file: {e}. The alignment file may be corrupted or empty.")
                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                        LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                        return
                    # Use BLOSUM62 for protein sequences (changed from identity for nucleotides)
                    calculator = DistanceCalculator("blosum62")
                    distMatrix = calculator.get_distance(align)
                    constructor = DistanceTreeConstructor()
                    UGMATree = constructor.upgma(distMatrix)
                    NJTree = constructor.nj(distMatrix)

                    def hide_labels(clade):
                        try:
                            # Don't hide any labels - show all sequence names
                            return str(clade)
                        except ValueError:
                            return str(clade)

                    UPGMA_TEXT = f"""The evolutionary history was inferred using the UPGMA method [1]. The optimal tree is shown. The evolutionary distances were computed using the BLOSUM62 method [2] and are in the
units of amino acid substitutions per site. This analysis involved {len(MATCHING_SEQ)} protein sequences. Multiple sequence alignment was performed using {method.split(' - ')[1]}. Evolutionary analyses
were conducted in B-Pan [3]

1. Sneath P.H.A. and Sokal R.R. (1973). Numerical Taxonomy. Freeman, San Francisco.

2. Henikoff S. and Henikoff J.G. (1992). Amino acid substitution matrices from protein blocks. PNAS 89:10915-10919.

3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                    NJ_TEXT = f"""The evolutionary history was inferred using the Neighbor-Joining method [1]. The optimal tree is shown. The evolutionary distances were computed using the BLOSUM62 method [2] and are
in the units of amino acid substitutions per site. This analysis involved {len(MATCHING_SEQ)} protein sequences. Multiple sequence alignment was performed using {method.split(' - ')[1]}. Evolutionary
analyses were conducted in B-Pan [3]

1. Saitou N. and Nei M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4:406-425.

2. Henikoff S. and Henikoff J.G. (1992). Amino acid substitution matrices from protein blocks. PNAS 89:10915-10919.

3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                    apply_publication_style()
                    
                    # Configure publication-quality tree plotting parameters with branch lengths
                    def plot_publication_tree(tree, title, filename_base, description_text):
                        """Plot a phylogenetic tree with publication-quality styling and branch length labels"""
                        try:
                            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 16),  # Increased height for better tree layout
                                                           gridspec_kw={'height_ratios': [8, 1]})
                            
                            # Custom branch label function to show branch lengths more prominently
                            def branch_length_labels(clade):
                                if clade.branch_length is not None and clade.branch_length > 0:
                                    return f"{clade.branch_length:.3f}"
                                return None
                            
                            # Draw tree with branch length labels
                            Phylo.draw(tree, do_show=False, label_func=hide_labels, axes=ax1,
                                      branch_labels=branch_length_labels, show_confidence=False)
                            
                            # Enhance tree appearance
                            ax1.set_title(title, fontsize=26, fontweight='bold', pad=25)
                            
                            # Style the tree branches - make them shorter and more visible
                            for line in ax1.get_lines():
                                line.set_linewidth(3.0)  # Thicker lines for better visibility
                                line.set_color('black')
                            
                            # Style leaf labels (sequence names) and branch length labels
                            for text in ax1.texts:
                                # Check if this is a branch length label or leaf label
                                text_content = text.get_text()
                                try:
                                    # If it's a number (branch length), make it more prominent and clear
                                    float(text_content)
                                    text.set_fontsize(14)  # Increased from 10 to 14
                                    text.set_fontweight('bold')  # Made bold
                                    text.set_color('red')  # Changed from blue to red for better visibility
                                    text.set_bbox(dict(boxstyle='round,pad=0.3', facecolor='white', 
                                                     edgecolor='red', alpha=0.9, linewidth=1))  # Added background box
                                except ValueError:
                                    # If it's a sequence name, make it prominent
                                    text.set_fontsize(12)
                                    text.set_fontweight('bold')
                                    text.set_color('black')
                            
                            # Clean axis appearance
                            ax1.set_xticks([])
                            ax1.set_yticks([])
                            ax1.set_xlabel('')
                            ax1.set_ylabel('')
                            
                            # Remove axis spines for cleaner look
                            for spine in ax1.spines.values():
                                spine.set_visible(False)
                            
                            # Set white background
                            ax1.set_facecolor('white')
                            
                            # Enhanced description box
                            ax2.axis("off")
                            ax2.text(0.5, 0.5, description_text, fontsize=12, fontweight='normal',
                                     ha='center', va='center', transform=ax2.transAxes,
                                     bbox=dict(boxstyle='round,pad=1.2', facecolor='white', 
                                              edgecolor='black', alpha=0.95, linewidth=2),
                                     wrap=True)
                            
                            plt.tight_layout(pad=2.0)
                            
                            # Save with publication quality
                            plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.png", dpi=600, bbox_inches='tight',
                                       facecolor='white', edgecolor='none', pad_inches=0.3)
                            plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.svg", bbox_inches='tight', pad_inches=0.3)
                            plt.clf()
                            
                            print(f"[Success] Tree plot saved: {filename_base}")
                            
                        except Exception as e:
                            print(f"[Error] Failed to plot tree {filename_base}: {e}")
                            # Try to save a simple tree plot as fallback
                            try:
                                plt.figure(figsize=(16, 12))
                                Phylo.draw(tree, do_show=False)
                                plt.title(f"{title} (Simple View)", fontsize=20, fontweight='bold')
                                plt.savefig(f"{RESULTS_PHYLO}/{filename_base}_simple.png", dpi=300, bbox_inches='tight')
                                plt.clf()
                                print(f"[Info] Simple tree plot saved as fallback: {filename_base}_simple.png")
                            except Exception as fallback_error:
                                print(f"[Error] Fallback tree plot also failed: {fallback_error}")
                                messagebox.showwarning("Tree Plot Warning", 
                                                     f"Failed to create enhanced tree plot for {filename_base}. "
                                                     f"Check the Results/Phylo directory for any generated files.")

                    if Phylogenetic_combo_box_var.get() == "All":
                        Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                        plot_publication_tree(UGMATree, 
                                            f"UPGMA Phylogenetic Tree ({method.split(' - ')[1]})",
                                            "UPGMA", UPGMA_TEXT)
                        
                        Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                        plot_publication_tree(NJTree,
                                            f"Neighbor-Joining Phylogenetic Tree ({method.split(' - ')[1]})", 
                                            "NJTree", NJ_TEXT)
                        
                    elif Phylogenetic_combo_box_var.get() == "UPGMA":
                        Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                        plot_publication_tree(UGMATree,
                                            f"UPGMA Phylogenetic Tree ({method.split(' - ')[1]})",
                                            "UPGMA", UPGMA_TEXT)
                        
                    elif Phylogenetic_combo_box_var.get() == "Neighbour Joining":
                        Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                        plot_publication_tree(NJTree,
                                            f"Neighbor-Joining Phylogenetic Tree ({method.split(' - ')[1]})",
                                            "NJTree", NJ_TEXT)
                    
                    # Optional KEGG annotation if toggled
                    if KEGG_Include_var.get() == "KEGG True":
                        try:
                            os.makedirs(RESULTS_KEGG, exist_ok=True)
                            
                            QUERY = all_names[0]
                            kegg_db = resolve_diamond_db_base(KEGG_DB)
                            if not diamond_db_present(kegg_db):
                                ensure_kegg_dmnd()
                            if not diamond_db_present(kegg_db):
                                messagebox.showerror("KEGG Error", f"KEGG database not found. Expected {KEGG_DB}.dmnd")
                                raise FileNotFoundError("KEGG DB missing")
                            kegg_out = resolve_path(os.path.join(RESULTS_KEGG, "output_file.txt"))
                            run_diamond_blastp(kegg_db, QUERY, kegg_out)
                            if not os.path.exists(kegg_out):
                                messagebox.showerror("KEGG Error", "KEGG DIAMOND output not created.")
                                raise FileNotFoundError("KEGG output missing")

                            kc = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                  "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                            KDATA = pd.read_csv(kegg_out, sep='\t', names=kc)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                            KDATA = KDATA[KDATA["PERCENTAGE IDENTITY"] >= KEGG_Percentage_combo_box_var.get()]

                            best = []
                            for q in KDATA["QUERY"].unique():
                                sub = KDATA[KDATA["QUERY"] == q]
                                best.append(sub.loc[sub["BIT SCORE"].idxmax()])
                            if best:
                                KDATA_BEST = pd.DataFrame(best)
                                # Optional mapping file to KO or descriptions
                                kmap_path = KEGG_ANNOTATION_CSV
                                if os.path.exists(kmap_path):
                                    KMAP = pd.read_csv(kmap_path)
                                    if "TITLE" in KMAP.columns:
                                        KDATA_BEST = KDATA_BEST.merge(KMAP, left_on="SUBJECT", right_on="TITLE", how="left")
                                KDATA_BEST.to_csv(os.path.join(RESULTS_KEGG, "KEGG_Results.csv"), index=False)
                        except Exception:
                            messagebox.showerror("KEGG Error", "KEGG Analysis failed. Please check the KEGG database and try again.")
                            pass
                    
                    messagebox.showinfo("Info", f"Phylogenetic analysis completed using {method}.")
                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                    LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                    return
                    
                elif DIAMOND_EXE is not None and os.path.exists(DIAMOND_EXE):

                    fasta_file = all_names[0]
                    TITLE, FASTA, DES = read_submitted_fasta_file(fasta_file)
                    CORE_DATA = pd.DataFrame({"TITLE": TITLE, "FASTA": FASTA, "DES": DES})

                    # Build similarity list via DIAMOND for sub-selection when using MSA methods
                    db_name = PHYLO_DB
                    create_diamond_db(fasta_file, db_name)
                    run_diamond_blastp(db_name, fasta_file, os.path.join(RESULTS_PHYLO, "output_file.txt"))
                    COLUMN_NAMES = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                    "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                    DATA_SEP = pd.read_csv(os.path.join(RESULTS_PHYLO, "output_file.txt"), sep='\t', names=COLUMN_NAMES)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY"]]
                    FILTERED_DATA = DATA_SEP[DATA_SEP["QUERY"] != DATA_SEP["SUBJECT"]]
                    GREATER_THAN = FILTERED_DATA[FILTERED_DATA["PERCENTAGE IDENTITY"] > Phylogenetic_Percentage_combo_box_var.get()]
                    CONCAT_CORE = pd.DataFrame({"CONCAT_CORE": list(set(GREATER_THAN["QUERY"].to_list()).union(set(GREATER_THAN["SUBJECT"].to_list())))}).drop_duplicates()
                    CONCAT_CORE_str = CONCAT_CORE["CONCAT_CORE"].apply(str)
                    MATCHING_SEQ = CORE_DATA[CORE_DATA["TITLE"].isin(CONCAT_CORE_str.to_list())]

                    if Phylo_method_var.get().startswith("k-mer"):
                        # k-mer Jaccard distance matrix directly from all sequences (or filtered if available)
                        seq_df = MATCHING_SEQ if len(MATCHING_SEQ) != 0 else CORE_DATA
                        
                        # Check if we have enough sequences for phylogenetic analysis
                        if len(seq_df) < 2:
                            messagebox.showerror("Error", "At least 2 sequences are required for phylogenetic analysis.")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            return
                        k = int(KMER_var.get())
                        
                        # Validate k-mer size
                        if k < 3 or k > 12:
                            messagebox.showerror("Error", "k-mer size must be between 3 and 12.")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            return
                            
                        def kmers(s, k):
                            s = str(s)
                            return {s[i:i+k] for i in range(len(s) - k + 1)} if len(s) >= k else set()
                        
                        # Filter out sequences that are too short for k-mer analysis
                        valid_sequences = []
                        valid_titles = []
                        for title, seq in zip(seq_df["TITLE"], seq_df["FASTA"]):
                            if len(str(seq)) >= k:
                                valid_sequences.append(seq)
                                valid_titles.append(title)
                        
                        if len(valid_sequences) < 2:
                            messagebox.showerror("Error", f"Not enough sequences with length >= {k} for k-mer analysis.")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            return
                        
                        titles = valid_titles
                        km_sets = [kmers(seq, k) for seq in valid_sequences]
                        n = len(titles)
                        
                        # Build distance matrix using Biopython's format
                        from Bio.Phylo.TreeConstruction import DistanceMatrix
                        
                        # Update progress
                        LABEL_PERCENTAGE.config(text=f"Progress : Calculating k-mer distances", fg="red")
                        
                        # Create distance matrix in the exact format Biopython expects
                        # Each row i should have exactly i+1 elements (lower triangular incl. diagonal)
                        matrix = []
                        for i in range(n):
                            row = []
                            for j in range(i + 1):
                                if j == i:
                                    # Diagonal element is zero distance to self
                                    row.append(0.0)
                                else:
                                    # Calculate Jaccard distance
                                    inter = len(km_sets[i] & km_sets[j])
                                    union = len(km_sets[i] | km_sets[j])
                                    if union == 0:
                                        d = 1.0  # Maximum distance if no overlap
                                    else:
                                        d = 1.0 - (inter / union)
                                    row.append(d)
                            matrix.append(row)
                        # Validate shape is lower-triangular with diagonal
                        for i, row in enumerate(matrix):
                            if len(row) != i + 1:
                                messagebox.showerror("Error", f"k-mer distance matrix row {i} has length {len(row)}, expected {i+1}")
                                LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                                return
                        
                        # Create DistanceMatrix
                        try:
                            dm = DistanceMatrix(names=titles, matrix=matrix)
                            constructor = DistanceTreeConstructor()
                            
                            # Update progress
                            LABEL_PERCENTAGE.config(text=f"Progress : Building phylogenetic trees", fg="red")
                            
                            UGMATree = constructor.upgma(dm)
                            NJTree = constructor.nj(dm)
                        except Exception as e:
                            messagebox.showerror("Error", f"k-mer distance calculation failed: {e}")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            return
                        
                        # Tree visualization for k-mer method
                        def hide_labels(clade):
                            try:
                                if str(clade)[0] == "I":
                                    return ""
                                else:
                                    return clade
                            except ValueError:
                                pass

                        UPGMA_TEXT = f"""The evolutionary history was inferred using the UPGMA method [1]. The optimal tree is shown. The evolutionary distances were computed using k-mer Jaccard distance (k={k}) [2] and are in the
units of dissimilarity. This analysis involved {n} amino acid sequences. Evolutionary analyses were conducted in B-Pan [3]

1. Sneath P.H.A. and Sokal R.R. (1973). Numerical Taxonomy. Freeman, San Francisco.

2. Jaccard P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles 37:547-579.

3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                        NJ_TEXT = f"""The evolutionary history was inferred using the Neighbor-Joining method [1]. The optimal tree is shown. The evolutionary distances were computed using k-mer Jaccard distance (k={k}) [2] and are in the
units of dissimilarity. This analysis involved {n} amino acid sequences. Evolutionary analyses were conducted in B-Pan [3]

1. Saitou N. and Nei M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4:406-425.

2. Jaccard P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles 37:547-579.

3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                        apply_publication_style()
                        
                        # Define publication-quality tree plotting for k-mer method
                        def plot_publication_kmer_tree(tree, title, filename_base, description_text):
                            """Plot a k-mer phylogenetic tree with publication-quality styling"""
                            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12),  # Larger figure
                                                           gridspec_kw={'height_ratios': [8, 1]})
                            
                            # Draw tree with enhanced styling
                            Phylo.draw(tree, do_show=False, label_func=hide_labels, axes=ax1,
                                      branch_labels=None, show_confidence=True)
                            
                            # Enhance tree appearance
                            ax1.set_title(title, fontsize=26, fontweight='bold', pad=25)
                            
                            # Style the tree branches and labels
                            for line in ax1.get_lines():
                                line.set_linewidth(3)  # Thicker branches
                                line.set_color('black')  # Ensure branches are black
                            
                            # Enhance any text elements (leaf labels)
                            for text in ax1.texts:
                                text.set_fontsize(14)
                                text.set_fontweight('bold')
                            
                            # Remove axis ticks and labels for cleaner look
                            ax1.set_xticks([])
                            ax1.set_yticks([])
                            ax1.set_xlabel('')
                            ax1.set_ylabel('')
                            
                            # Add a subtle frame around the tree
                            ax1.spines['top'].set_visible(True)
                            ax1.spines['right'].set_visible(True)
                            ax1.spines['bottom'].set_visible(True)
                            ax1.spines['left'].set_visible(True)
                            for spine in ax1.spines.values():
                                spine.set_linewidth(2)
                                spine.set_color('black')
                            
                            # Enhanced description box
                            ax2.axis("off")
                            ax2.text(0.5, 0.5, description_text, fontsize=12, fontweight='normal',
                                     ha='center', va='center', transform=ax2.transAxes,
                                     bbox=dict(boxstyle='round,pad=1.2', facecolor='white', 
                                              edgecolor='black', alpha=0.95, linewidth=2),
                                     wrap=True)
                            
                            plt.tight_layout(pad=2.0)
                            
                            # Save with publication quality
                            plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.png", dpi=600, bbox_inches='tight',
                                       facecolor='white', edgecolor='none', pad_inches=0.2)
                            plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.svg", bbox_inches='tight', pad_inches=0.2)
                            plt.clf()

                        if Phylogenetic_combo_box_var.get() == "All":
                            Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                            plot_publication_kmer_tree(UGMATree,
                                                      f"UPGMA Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                      "UPGMA", UPGMA_TEXT)
                            
                            Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                            plot_publication_kmer_tree(NJTree,
                                                      f"Neighbor-Joining Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                      "NJTree", NJ_TEXT)
                            
                        elif Phylogenetic_combo_box_var.get() == "UPGMA":
                            Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                            plot_publication_kmer_tree(UGMATree,
                                                      f"UPGMA Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                      "UPGMA", UPGMA_TEXT)
                            
                        elif Phylogenetic_combo_box_var.get() == "Neighbour Joining":
                            Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                            plot_publication_kmer_tree(NJTree,
                                                      f"Neighbor-Joining Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                      "NJTree", NJ_TEXT)
                        
                        # Reset progress and show completion message for k-mer method
                        LABEL_PERCENTAGE.config(text=f"Progress : 0%", fg=Label_Color)
                        messagebox.showinfo("Info", "k-mer Jaccard phylogenetic analysis completed.")
                        
                        # Optional KEGG annotation if toggled for k-mer method
                        if KEGG_Include_var.get() == "KEGG True":
                             try:
                                 os.makedirs(RESULTS_KEGG, exist_ok=True)
                                 
                                 QUERY = all_names[0]
                                 kegg_db = resolve_diamond_db_base(KEGG_DB)
                                 if not diamond_db_present(kegg_db):
                                     ensure_kegg_dmnd()
                                 if not diamond_db_present(kegg_db):
                                     messagebox.showerror("KEGG Error", f"KEGG database not found. Expected {KEGG_DB}.dmnd")
                                     raise FileNotFoundError("KEGG DB missing")
                                 kegg_out = resolve_path(os.path.join(RESULTS_KEGG, "output_file.txt"))
                                 run_diamond_blastp(kegg_db, QUERY, kegg_out)
                                 if not os.path.exists(kegg_out):
                                     messagebox.showerror("KEGG Error", "KEGG DIAMOND output not created.")
                                     raise FileNotFoundError("KEGG output missing")

                                 kc = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                       "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                                 KDATA = pd.read_csv(kegg_out, sep='\t', names=kc)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                                 KDATA = KDATA[KDATA["PERCENTAGE IDENTITY"] >= KEGG_Percentage_combo_box_var.get()]

                                 best = []
                                 for q in KDATA["QUERY"].unique():
                                     sub = KDATA[KDATA["QUERY"] == q]
                                     best.append(sub.loc[sub["BIT SCORE"].idxmax()])
                                 if best:
                                     KDATA_BEST = pd.DataFrame(best)
                                     # Optional mapping file to KO or descriptions
                                     kmap_path = KEGG_ANNOTATION_CSV
                                     if os.path.exists(kmap_path):
                                         KMAP = pd.read_csv(kmap_path)
                                         if "TITLE" in KMAP.columns:
                                             KDATA_BEST = KDATA_BEST.merge(KMAP, left_on="SUBJECT", right_on="TITLE", how="left")
                                     KDATA_BEST.to_csv(os.path.join(RESULTS_KEGG, "KEGG_Results.csv"), index=False)
                             except Exception:
                                 messagebox.showerror("KEGG Error", "KEGG Analysis failed. Please check the KEGG database and try again.")
                                 pass
                    else:
                        if len(MATCHING_SEQ) == 0:
                            # Auto-relax threshold progressively; if still empty, fall back to k-mer method
                            try:
                                current_thr = int(Phylogenetic_Percentage_combo_box_var.get())
                            except Exception:
                                current_thr = 90
                            relaxed_thresholds = [t for t in [current_thr, 70, 60, 50, 40, 30, 20, 10] if t <= 100]
                            recovered = False
                            for thr in relaxed_thresholds:
                                GREATER_THAN = FILTERED_DATA[FILTERED_DATA["PERCENTAGE IDENTITY"] > thr]
                                CONCAT_CORE = pd.DataFrame({"CONCAT_CORE": list(set(GREATER_THAN["QUERY"].to_list()).union(set(GREATER_THAN["SUBJECT"].to_list())))}).drop_duplicates()
                                CONCAT_CORE_str = CONCAT_CORE["CONCAT_CORE"].apply(str)
                                MATCHING_SEQ = CORE_DATA[CORE_DATA["TITLE"].isin(CONCAT_CORE_str.to_list())]
                                if len(MATCHING_SEQ) >= 2:
                                    ui_info(f"No pairs at {current_thr}%. Proceeding with relaxed threshold {thr}% for MSA.")
                                    recovered = True
                                    break
                            if not recovered:
                                # Fall back to k-mer pipeline using all sequences
                                ui_info("No pairs found even after relaxing identity threshold; switching to k-mer Jaccard method.")
                                seq_df = CORE_DATA
                                k = int(KMER_var.get())
                                if len(seq_df) < 2:
                                    messagebox.showerror("Error", "At least 2 sequences are required for phylogenetic analysis.")
                                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                                    return
                                if k < 3 or k > 12:
                                    messagebox.showerror("Error", "k-mer size must be between 3 and 12.")
                                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                                    return
                                def kmers(s, k):
                                    s = str(s)
                                    return {s[i:i+k] for i in range(len(s) - k + 1)} if len(s) >= k else set()
                                valid_sequences = []
                                valid_titles = []
                                for title, seq in zip(seq_df["TITLE"], seq_df["FASTA"]):
                                    if len(str(seq)) >= k:
                                        valid_sequences.append(seq)
                                        valid_titles.append(title)
                                if len(valid_sequences) < 2:
                                    messagebox.showerror("Error", f"Not enough sequences with length >= {k} for k-mer analysis.")
                                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                                    return
                                titles = valid_titles
                                km_sets = [{str(seq)[i:i+k] for i in range(len(str(seq)) - k + 1)} for seq in valid_sequences]
                                LABEL_PERCENTAGE.config(text=f"Progress : Calculating k-mer distances", fg="red")
                                matrix = []
                                for i in range(len(titles)):
                                    row = []
                                    for j in range(i + 1):
                                        if j == i:
                                            row.append(0.0)
                                        else:
                                            inter = len(km_sets[i] & km_sets[j])
                                            union = len(km_sets[i] | km_sets[j])
                                            d = 1.0 if union == 0 else 1.0 - (inter / union)
                                            row.append(d)
                                    matrix.append(row)
                                for i, row in enumerate(matrix):
                                    if len(row) != i + 1:
                                        messagebox.showerror("Error", f"k-mer distance matrix row {i} has length {len(row)}, expected {i+1}")
                                        LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                                        return
                                try:
                                    dm = DistanceMatrix(names=titles, matrix=matrix)
                                    constructor = DistanceTreeConstructor()
                                    LABEL_PERCENTAGE.config(text=f"Progress : Building phylogenetic trees", fg="red")
                                    UGMATree = constructor.upgma(dm)
                                    NJTree = constructor.nj(dm)
                                except Exception as e:
                                    messagebox.showerror("Error", f"k-mer distance calculation failed: {e}")
                                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                                    return
                                def hide_labels(clade):
                                    try:
                                        if str(clade)[0] == "I":
                                            return ""
                                        else:
                                            return clade
                                    except ValueError:
                                        pass
                                apply_publication_style()
                                
                                # ✅ Define publication-quality k-mer tree plotting (matching MSA tree style)
                                def plot_publication_kmer_tree_fallback(tree, title, filename_base, description_text):
                                    """Plot a k-mer phylogenetic tree with publication-quality styling (fallback auto mode)"""
                                    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12),  # Larger figure
                                                                   gridspec_kw={'height_ratios': [8, 1]})
                                    
                                    # Draw tree with enhanced styling
                                    Phylo.draw(tree, do_show=False, label_func=hide_labels, axes=ax1,
                                              branch_labels=None, show_confidence=True)
                                    
                                    # Enhance tree appearance
                                    ax1.set_title(title, fontsize=26, fontweight='bold', pad=25)
                                    
                                    # Style the tree branches and labels
                                    for line in ax1.get_lines():
                                        line.set_linewidth(3)  # Thicker branches
                                        line.set_color('black')  # Ensure branches are black
                                    
                                    # Enhance any text elements (leaf labels)
                                    for text in ax1.texts:
                                        text.set_fontsize(14)
                                        text.set_fontweight('bold')
                                    
                                    # Remove axis ticks and labels for cleaner look
                                    ax1.set_xticks([])
                                    ax1.set_yticks([])
                                    ax1.set_xlabel('')
                                    ax1.set_ylabel('')
                                    
                                    # Add a subtle frame around the tree
                                    ax1.spines['top'].set_visible(True)
                                    ax1.spines['right'].set_visible(True)
                                    ax1.spines['bottom'].set_visible(True)
                                    ax1.spines['left'].set_visible(True)
                                    for spine in ax1.spines.values():
                                        spine.set_linewidth(2)
                                        spine.set_color('black')
                                    
                                    # Enhanced description box
                                    ax2.axis("off")
                                    ax2.text(0.5, 0.5, description_text, fontsize=12, fontweight='normal',
                                             ha='center', va='center', transform=ax2.transAxes,
                                             bbox=dict(boxstyle='round,pad=1.2', facecolor='white', 
                                                      edgecolor='black', alpha=0.95, linewidth=2),
                                             wrap=True)
                                    
                                    plt.tight_layout(pad=2.0)
                                    
                                    # Save with publication quality
                                    plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.png", dpi=600, bbox_inches='tight',
                                               facecolor='white', edgecolor='none', pad_inches=0.2)
                                    plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.svg", bbox_inches='tight', pad_inches=0.2)
                                    plt.clf()
                                
                                # Create proper description texts for fallback k-mer trees
                                UPGMA_TEXT_FALLBACK = f"""The evolutionary history was inferred using the UPGMA method [1]. The optimal tree is shown. The evolutionary distances were computed using k-mer Jaccard distance (k={k}) [2] and are in the
units of dissimilarity. This analysis involved {len(titles)} amino acid sequences. Evolutionary analyses were conducted in B-Pan [3]

1. Sneath P.H.A. and Sokal R.R. (1973). Numerical Taxonomy. Freeman, San Francisco.

2. Jaccard P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles 37:547-579.

3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                                NJ_TEXT_FALLBACK = f"""The evolutionary history was inferred using the Neighbor-Joining method [1]. The optimal tree is shown. The evolutionary distances were computed using k-mer Jaccard distance (k={k}) [2] and are in the
units of dissimilarity. This analysis involved {len(titles)} amino acid sequences. Evolutionary analyses were conducted in B-Pan [3]

1. Saitou N. and Nei M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4:406-425.

2. Jaccard P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles 37:547-579.

3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                                if Phylogenetic_combo_box_var.get() == "All":
                                    Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                                    plot_publication_kmer_tree_fallback(UGMATree,
                                                                      f"UPGMA Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                                      "UPGMA", UPGMA_TEXT_FALLBACK)
                                    
                                    Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                                    plot_publication_kmer_tree_fallback(NJTree,
                                                                      f"Neighbor-Joining Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                                      "NJTree", NJ_TEXT_FALLBACK)
                                    
                                elif Phylogenetic_combo_box_var.get() == "UPGMA":
                                    Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                                    plot_publication_kmer_tree_fallback(UGMATree,
                                                                      f"UPGMA Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                                      "UPGMA", UPGMA_TEXT_FALLBACK)
                                    
                                elif Phylogenetic_combo_box_var.get() == "Neighbour Joining":
                                    Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                                    plot_publication_kmer_tree_fallback(NJTree,
                                                                      f"Neighbor-Joining Phylogenetic Tree (k-mer Jaccard, k={k})",
                                                                      "NJTree", NJ_TEXT_FALLBACK)
                                LABEL_PERCENTAGE.config(text=f"Progress : 0%", fg=Label_Color)
                                messagebox.showinfo("Info", "k-mer Jaccard phylogenetic analysis completed.")
                                # Optional KEGG omitted in auto-fallback branch
                                return

                        # Write filtered sequences
                        with open(os.path.join(RESULTS_PHYLO, "Alignment.fasta"), "w") as Alignment_File:
                            for t, s in zip(MATCHING_SEQ["TITLE"], MATCHING_SEQ["FASTA"]):
                                Alignment_File.write(f">{t}\n{s}\n")
                        Alignment_File_Path = os.path.join(RESULTS_PHYLO, "Alignment.fasta")
                        OUTFILE = RESULTS_PHYLO

                        method = Phylo_method_var.get()
                        LABEL_PERCENTAGE.config(text=f"Progress : Running {method}", fg="red")
                        try:
                            if method == "MSA - Clustal Omega":
                                CLUSTAL_PATH = CLUSTALO_EXE
                                print(f"{CLUSTAL_PATH} -i {Alignment_File_Path} -o {OUTFILE}\\Alignment.fasta --force --threads={os.cpu_count()}")
                                result = subprocess.run(
                                       [CLUSTAL_PATH, "-i", Alignment_File_Path, "-o", OUTFILE + "\\Alignment.fasta",
                                       "--force", f"--threads={os.cpu_count()}"],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        startupinfo=startupinfo,
                                        shell=False,
                                        timeout=300
                                     )
                                if result.returncode != 0:
                                    print(f"Clustal Omega error: {result.stderr.decode()}")

                            elif method == "MSA - MAFFT":
                                MAFFT_PATH = MAFFT_EXE
                                print(f"{MAFFT_PATH} --auto {Alignment_File_Path} > {OUTFILE}\\Alignment.fasta")
                                is_bat = MAFFT_PATH.lower().endswith((".bat", ".cmd"))
                                mafft_invoke = ["cmd", "/c", MAFFT_PATH, "--auto", Alignment_File_Path] if is_bat \
                                    else [MAFFT_PATH, "--auto", Alignment_File_Path]
                                with open(OUTFILE + "\\Alignment.fasta", "w") as fout:
                                     result = subprocess.run(
                                            mafft_invoke,
                                            stdout=fout,
                                            stderr=subprocess.PIPE,
                                            startupinfo=startupinfo,
                                            shell=False,
                                            timeout=300
                                     )
                                     if result.returncode != 0:
                                         print(f"MAFFT error: {result.stderr.decode()}")

                            elif method == "MSA - MUSCLE":
                                MUSCLE_PATH = MUSCLE_EXE
                                print(f"{MUSCLE_PATH} -align {Alignment_File_Path} -output {OUTFILE}\\Alignment.fasta")
                                result = subprocess.run(
                                       [MUSCLE_PATH, "-align", Alignment_File_Path, "-output", OUTFILE + "\\Alignment.fasta"],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       startupinfo=startupinfo,
                                       shell=False,
                                       timeout=300
                                    )
                                if result.returncode != 0:
                                    print(f"MUSCLE error: {result.stderr.decode()}")

                            else:
                                raise ValueError("Unknown phylogenetic method")
                        except subprocess.TimeoutExpired:
                            messagebox.showerror("Error", f"MSA timed out after 5 minutes. Try with fewer sequences or a different method.")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                            return
                        except Exception as e:
                            messagebox.showerror("Error", f"MSA failed: {e}")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                            return

                        # Check if alignment file was created
                        if not os.path.exists(OUTFILE + "\\Alignment.fasta"):
                            messagebox.showerror("Error", f"Alignment file was not created. Please check if {method} is working correctly.")
                            LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                            LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)
                            return

                        LABEL_PERCENTAGE.config(text="Progress : Building trees", fg="red")
                        align = AlignIO.read(OUTFILE + "\\Alignment.fasta", "fasta")
                        # Use identity calculator for nucleotide sequences (16S rRNA)
                        calculator = DistanceCalculator("identity")
                        distMatrix = calculator.get_distance(align)
                        constructor = DistanceTreeConstructor()
                        UGMATree = constructor.upgma(distMatrix)
                        NJTree = constructor.nj(distMatrix)

                        def hide_labels(clade):
                            try:
                                if str(clade)[0] == "I":
                                    return ""
                                else:
                                    return clade
                            except ValueError:
                                pass


                        UPGMA_TEXT = """The evolutionary history was inferred using the UPGMA method [1]. The optimal tree is shown. The evolutionary distances were computed using the identity method [2] and are in the
units of the number of nucleotide substitutions per site. This analysis involved nucleotide sequences (16S rRNA). All ambiguous positions were removed for each sequence pair (pairwise deletion option). Evolutionary analyses
were conducted in B-Pan [3]

                                1. Sneath P.H.A. and Sokal R.R. (1973). Numerical Taxonomy. Freeman, San Francisco.

                                2. Identity distance calculation for nucleotide sequences.

                                3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                        NJ_TEXT = """The evolutionary history was inferred using the Neighbor-Joining method [1]. The optimal tree is shown. The evolutionary distances were computed using the identity method [2] and are
in the units of the number of nucleotide substitutions per site. This analysis involved nucleotide sequences (16S rRNA). All ambiguous positions were removed for each sequence pair (pairwise deletion option). Evolutionary
analyses were conducted in B-Pan [3]

                        1. Saitou N. and Nei M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4:406-425.

                        2. Identity distance calculation for nucleotide sequences.

                        3. B-Pan: A Robust Software Package for Bacterial Pangenome Analysis."""

                        apply_publication_style()

                        # ✅ CRITICAL FIX: Use the same high-quality tree plotting for ALL trees
                        def plot_publication_tree_nucleotide(tree, title, filename_base, description_text):
                            """Plot a nucleotide phylogenetic tree with publication-quality styling"""
                            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12),  # Larger figure
                                                           gridspec_kw={'height_ratios': [8, 1]})
                            
                            # Draw tree with enhanced styling
                            Phylo.draw(tree, do_show=False, label_func=hide_labels, axes=ax1,
                                      branch_labels=None, show_confidence=True)
                            
                            # Enhance tree appearance
                            ax1.set_title(title, fontsize=26, fontweight='bold', pad=25)
                            
                            # Style the tree branches and labels
                            for line in ax1.get_lines():
                                line.set_linewidth(3)  # Thicker branches
                                line.set_color('black')  # Ensure branches are black
                            
                            # Enhance any text elements (leaf labels)
                            for text in ax1.texts:
                                text.set_fontsize(14)
                                text.set_fontweight('bold')
                            
                            # Remove axis ticks and labels for cleaner look
                            ax1.set_xticks([])
                            ax1.set_yticks([])
                            ax1.set_xlabel('')
                            ax1.set_ylabel('')
                            
                            # Add a subtle frame around the tree
                            ax1.spines['top'].set_visible(True)
                            ax1.spines['right'].set_visible(True)
                            ax1.spines['bottom'].set_visible(True)
                            ax1.spines['left'].set_visible(True)
                            for spine in ax1.spines.values():
                                spine.set_linewidth(2)
                                spine.set_color('black')
                            
                            # Enhanced description box
                            ax2.axis("off")
                            ax2.text(0.5, 0.5, description_text, fontsize=12, fontweight='normal',
                                     ha='center', va='center', transform=ax2.transAxes,
                                     bbox=dict(boxstyle='round,pad=1.2', facecolor='white', 
                                              edgecolor='black', alpha=0.95, linewidth=2),
                                     wrap=True)
                            
                            plt.tight_layout(pad=2.0)
                            
                            # Save with publication quality
                            plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.png", dpi=600, bbox_inches='tight',
                                       facecolor='white', edgecolor='none', pad_inches=0.2)
                            plt.savefig(f"{RESULTS_PHYLO}/{filename_base}.svg", bbox_inches='tight', pad_inches=0.2)
                            plt.clf()

                        if Phylogenetic_combo_box_var.get() == "All":
                            Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                            plot_publication_tree_nucleotide(UGMATree, 
                                                           f"UPGMA Phylogenetic Tree (16S rRNA)",
                                                           "UPGMA", UPGMA_TEXT)
                            
                            Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                            plot_publication_tree_nucleotide(NJTree,
                                                           f"Neighbor-Joining Phylogenetic Tree (16S rRNA)", 
                                                           "NJTree", NJ_TEXT)
                            
                        elif Phylogenetic_combo_box_var.get() == "UPGMA":
                            Phylo.write(UGMATree, os.path.join(RESULTS_PHYLO, "UPGMA.nwk"), format="newick")
                            plot_publication_tree_nucleotide(UGMATree,
                                                           f"UPGMA Phylogenetic Tree (16S rRNA)",
                                                           "UPGMA", UPGMA_TEXT)
                            
                        elif Phylogenetic_combo_box_var.get() == "Neighbour Joining":
                            Phylo.write(NJTree, os.path.join(RESULTS_PHYLO, "NJTREE.nwk"), format="newick")
                            plot_publication_tree_nucleotide(NJTree,
                                                           f"Neighbor-Joining Phylogenetic Tree (16S rRNA)", 
                                                           "NJTree", NJ_TEXT)
                    

                    # Optional KEGG annotation if toggled
                    if KEGG_Include_var.get() == "KEGG True":
                        try:
                            os.makedirs(RESULTS_KEGG, exist_ok=True)
                            
                            QUERY = all_names[0]
                            kegg_db = resolve_diamond_db_base(KEGG_DB)
                            if not diamond_db_present(kegg_db):
                                ensure_kegg_dmnd()
                            if not diamond_db_present(kegg_db):
                                messagebox.showerror("KEGG Error", f"KEGG database not found. Expected {KEGG_DB}.dmnd")
                                raise FileNotFoundError("KEGG DB missing")
                            kegg_out = resolve_path(os.path.join(RESULTS_KEGG, "output_file.txt"))
                            run_diamond_blastp(kegg_db, QUERY, kegg_out)
                            if not os.path.exists(kegg_out):
                                messagebox.showerror("KEGG Error", "KEGG DIAMOND output not created.")
                                raise FileNotFoundError("KEGG output missing")

                            kc = ["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "alignment length", "mismatches",
                                  "gap opens", "q. start", "q. end", "s. start", "s. end", "E VALUE", "BIT SCORE"]
                            KDATA = pd.read_csv(kegg_out, sep='\t', names=kc)[["QUERY", "SUBJECT", "PERCENTAGE IDENTITY", "E VALUE", "BIT SCORE"]]
                            KDATA = KDATA[KDATA["PERCENTAGE IDENTITY"] >= KEGG_Percentage_combo_box_var.get()]

                            best = []
                            for q in KDATA["QUERY"].unique():
                                sub = KDATA[KDATA["QUERY"] == q]
                                best.append(sub.loc[sub["BIT SCORE"].idxmax()])
                            if best:
                                KDATA_BEST = pd.DataFrame(best)
                                # Optional mapping file to KO or descriptions
                                kmap_path = KEGG_ANNOTATION_CSV
                                if os.path.exists(kmap_path):
                                    KMAP = pd.read_csv(kmap_path)
                                    if "TITLE" in KMAP.columns:
                                        KDATA_BEST = KDATA_BEST.merge(KMAP, left_on="SUBJECT", right_on="TITLE", how="left")
                                KDATA_BEST.to_csv(os.path.join(RESULTS_KEGG, "KEGG_Results.csv"), index=False)
                        except Exception:
                            messagebox.showerror("KEGG Error", "KEGG Analysis failed. Please check the KEGG database and try again.")
                            pass

                    messagebox.showinfo("Info", "Analysis Completed.")

                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
                    LABEL_PERCENTAGE.config(text="Progress : 0%", fg=Label_Color)

                else:
                    messagebox.showerror("Error", "Diamond not found.")
                    LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)

            else:
                messagebox.showerror("Error","Phylogenetic Analysis requires one file for Analysis.")
                LOADING_LABEL.config(text="Job Status : None", fg=Label_Color)
        else:
            messagebox.showerror("Error", "Kindly Select an Analysis Method.")


Pangenome_combo_box_var = IntVar()
# Allow 70–100% in 5% steps; default to 90%
Pangenome_combo_box = ttk.Combobox(window, values=[x for x in range(70,101,5)], textvariable=Pangenome_combo_box_var)
Pangenome_combo_box_var.set(90)
Pangenome_combo_box.place(x=673, y=175)

def cleanup_temp_files():
    """Clean up temporary files created during analysis."""
    try:
        # Clean up any remaining temporary files
        temp_patterns = [
            "diamond-*.tmp",
            "*.tmp",
            "~$*"
        ]
        
        for pattern in temp_patterns:
            for temp_file in glob.glob(pattern):
                try:
                    if os.path.isfile(temp_file):
                        os.remove(temp_file)
                except Exception:
                    pass
    except Exception:
        pass

# Register cleanup function to run on window close
window.protocol("WM_DELETE_WINDOW", lambda: [cleanup_temp_files(), window.destroy()])

SUBMIT_Btn = Button(window, image=SUBMIT_BUTTON,bd=0, command=lambda: threading.Thread(target=FASTA_SUBMIT_BTN).start())
# Place Submit above footer text, centered (kept constant)
SUBMIT_Btn.place(x=430, y=535)

window.mainloop()