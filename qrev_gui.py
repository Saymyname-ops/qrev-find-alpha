"""
QREV GUI - Graphical Interface for Finding Best Alpha
Based on qrev_find_alpha.py but with tkinter GUI
"""

import numpy as np
import numba
from numba import njit
import sys
import csv
import time
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import queue

# Import the core functions from the original file
from qrev_find_alpha import (
    _qrev_sim, _qrev_rmse, load_experimental_data, 
    find_best_alpha, save_results, plot_results
)

class QREVGui:
    def __init__(self, root):
        self.root = root
        self.root.title("QREV - Find Best Alpha (GUI Version)")
        self.root.geometry("600x500")
        
        # Variables
        self.data_file = tk.StringVar()
        self.status_text = tk.StringVar(value="Ready")
        self.progress_var = tk.DoubleVar()
        self.is_running = False
        
        self.setup_ui()
        
    def setup_ui(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        
        # Title
        title_label = ttk.Label(main_frame, text="QREV - Find Best Alpha", 
                               font=("Arial", 16, "bold"))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # File selection
        ttk.Label(main_frame, text="Data File:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Entry(main_frame, textvariable=self.data_file, width=50).grid(
            row=1, column=1, sticky=(tk.W, tk.E), padx=(5, 5), pady=5)
        ttk.Button(main_frame, text="Browse...", 
                  command=self.browse_file).grid(row=1, column=2, padx=(5, 0), pady=5)
        
        # Parameters frame
        params_frame = ttk.LabelFrame(main_frame, text="Parameters", padding="10")
        params_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 10))
        params_frame.columnconfigure(1, weight=1)
        
        # Fixed parameters (matching original)
        ttk.Label(params_frame, text="Initial Potential:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Label(params_frame, text="0.0 V").grid(row=0, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(params_frame, text="Limit Potential:").grid(row=1, column=0, sticky=tk.W, pady=2)
        ttk.Label(params_frame, text="-30.0 V").grid(row=1, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(params_frame, text="Potential Step:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ttk.Label(params_frame, text="0.01 V").grid(row=2, column=1, sticky=tk.W, pady=2)
        
        # Control buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=10)
        
        self.run_button = ttk.Button(button_frame, text="Run Analysis", 
                                    command=self.run_analysis)
        self.run_button.pack(side=tk.LEFT, padx=(0, 10))
        
        self.stop_button = ttk.Button(button_frame, text="Stop", 
                                     command=self.stop_analysis, state="disabled")
        self.stop_button.pack(side=tk.LEFT, padx=(0, 10))
        
        ttk.Button(button_frame, text="Clear Results", 
                  command=self.clear_results).pack(side=tk.LEFT)
        
        # Progress frame
        progress_frame = ttk.LabelFrame(main_frame, text="Progress", padding="10")
        progress_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 10))
        progress_frame.columnconfigure(0, weight=1)
        
        ttk.Label(progress_frame, textvariable=self.status_text).grid(
            row=0, column=0, sticky=tk.W, pady=(0, 5))
        
        self.progress_bar = ttk.Progressbar(progress_frame, variable=self.progress_var, 
                                           maximum=100)
        self.progress_bar.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 5))
        
        self.progress_label = ttk.Label(progress_frame, text="0%")
        self.progress_label.grid(row=1, column=1, padx=(5, 0))
        
        # Results area
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="10")
        results_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 0))
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
        main_frame.rowconfigure(5, weight=1)
        
        # Text widget with scrollbar
        text_frame = ttk.Frame(results_frame)
        text_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        
        self.results_text = tk.Text(text_frame, height=15, width=70, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL, command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scrollbar.set)
        
        self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        # Queue for thread communication
        self.queue = queue.Queue()
        self.root.after(100, self.process_queue)
        
    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Choose Experimental Data File (Tab-separated: POT Y Q)",
            filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")],
            initialdir=r"C:\Users\ranin"
        )
        if filename:
            self.data_file.set(filename)
            
    def log_message(self, message):
        self.queue.put(("log", message))
        
    def update_status(self, message):
        self.queue.put(("status", message))
        
    def update_progress(self, value, text=None):
        self.queue.put(("progress", (value, text)))
        
    def process_queue(self):
        try:
            while True:
                msg_type, msg_data = self.queue.get_nowait()
                if msg_type == "log":
                    self.results_text.insert(tk.END, msg_data + "\n")
                    self.results_text.see(tk.END)
                elif msg_type == "status":
                    self.status_text.set(msg_data)
                elif msg_type == "progress":
                    value, text = msg_data
                    self.progress_var.set(value)
                    if text:
                        self.progress_label.config(text=text)
                    else:
                        self.progress_label.config(text=f"{int(value)}%")
        except queue.Empty:
            pass
        self.root.after(100, self.process_queue)
        
    def run_analysis(self):
        if self.is_running:
            return
            
        if not self.data_file.get():
            messagebox.showerror("Error", "Please select a data file first.")
            return
            
        # Disable UI during run
        self.is_running = True
        self.run_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.clear_results()
        
        # Start analysis in separate thread
        thread = threading.Thread(target=self.analysis_thread)
        thread.daemon = True
        thread.start()
        
    def stop_analysis(self):
        self.is_running = False
        self.update_status("Stopping...")
        self.stop_button.config(state="disabled")
        
    def clear_results(self):
        self.results_text.delete(1.0, tk.END)
        self.progress_var.set(0)
        self.progress_label.config(text="0%")
        self.status_text.set("Ready")
        
    def analysis_thread(self):
        try:
            self.log_message("=" * 60)
            self.log_message("QREV - FIND BEST ALPHA (GUI VERSION)")
            self.log_message("=" * 60)
            self.log_message(f"Data file: {self.data_file.get()}")
            
            # Fixed parameters (same as original)
            INIT = 0.0
            LIMIT = -30.0
            DELTA = 0.01
            
            self.log_message(f"INIT = {INIT}")
            self.log_message(f"LIMIT = {LIMIT}")
            self.log_message(f"DELTA = {DELTA}")
            
            # Alpha range (0.01 to 0.99)
            alphas = [round(i * 0.01, 4) for i in range(1, 100)]
            self.log_message(f"Alpha range: 0.01 to 0.99 ({len(alphas)} values)")
            
            # Omega values (focused on sensitive range)
            omegas = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 10**-1.5, 1e-1, 1e0, 1e1, 1e2, 1e3]
            self.log_message(f"Omega values: {len(omegas)}")
            self.log_message(f"Total combinations: {len(alphas) * len(omegas)}")
            self.log_message("=" * 60)
            
            # Load experimental data
            self.update_status("Loading experimental data...")
            experimental_data = load_experimental_data(self.data_file.get())
            if not experimental_data:
                self.log_message("ERROR: No data loaded from file.")
                return
                
            self.log_message(f"Loaded {len(experimental_data)} data points")
            self.log_message(f"  POT range: {experimental_data[0]['pot']:.2f} to {experimental_data[-1]['pot']:.2f}")
            y_values = [d['y'] for d in experimental_data]
            self.log_message(f"  Y range:   {min(y_values):.4e} to {max(y_values):.4e}")
            
            # Search
            self.log_message("Starting search...")
            all_results, best = find_best_alpha(
                experimental_data, alphas, omegas, INIT, LIMIT, DELTA
            )
            
            # Check if we should stop
            if not self.is_running:
                self.log_message("Analysis stopped by user.")
                return
                
            # Results
            self.log_message("=" * 60)
            self.log_message("BEST FIT FOUND")
            self.log_message("=" * 60)
            self.log_message(f"  Alpha:  {best['alpha']}")
            self.log_message(f"  Beta:   {best['beta']}")
            self.log_message(f"  Omega:  {best['omega']}")
            self.log_message(f"  RMSE:   {best['rmse']:.6e}")
            self.log_message("=" * 60)
            
            # Save results
            self.update_status("Saving results...")
            output_dir = Path(r"C:\Users\ranin\qrev_output_gui")
            save_results(all_results, best, experimental_data, INIT, LIMIT, DELTA, output_dir)
            
            # Plot results
            self.update_status("Generating plot...")
            plot_results(best, experimental_data, INIT, LIMIT, DELTA, output_dir)
            
            self.log_message("Analysis completed successfully!")
            self.update_status("Completed")
            
        except Exception as e:
            self.log_message(f"ERROR: {str(e)}")
            self.update_status("Error occurred")
        finally:
            # Re-enable UI
            self.is_running = False
            self.run_button.config(state="normal")
            self.stop_button.config(state="disabled")
            self.progress_var.set(0)
            self.progress_label.config(text="0%")

def main():
    root = tk.Tk()
    app = QREVGui(root)
    root.mainloop()

if __name__ == "__main__":
    main()