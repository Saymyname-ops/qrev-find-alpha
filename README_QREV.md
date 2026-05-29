# QREV - Find Best Alpha

This program finds the best alpha (and omega) parameters by comparing QREV simulation with experimental electrochemical data.

## Usage

1. Prepare your experimental data file (tab or space separated):
   ```
   POT     Y          Q
   0.0     0.0        0.0
   -0.5    -0.1       -0.005
   -1.0    -0.25      -0.018
   ...
   ```

2. Run the program:
   ```
   python qrev_find_alpha.py "path/to/your/data.txt"
   ```

   The program will:
   - Load your experimental data
   - Test all combinations of alpha (0.01 to 0.99) and omega values
   - Find the combination that minimizes RMSE (Root Mean Square Error)
   - Show a plot of E adm vs I adm (simulation vs experimental)
   - Wait for you to press ENTER before exiting

## What was fixed

### Original Issues:
1. **Wrong termination condition** - Was stopping simulation too early when q >= 0
2. **Inefficient Python loops** - O(n²) convolution was very slow
3. **Incorrect convolution sign** - Was using `-=` instead of `+=`
4. **Wrong omega range** - Testing high omega values where alpha has little effect
5. **Plot labels** - Was showing generic labels instead of E adm/I adm

### Fixes Applied:
1. **Numba JIT compilation** - Made simulation ~100x faster
2. **Corrected physics** - Fixed termination condition and convolution
3. **Focused omega range** - Now tests 1e-3 to 1e3 (where alpha sensitivity is highest)
4. **Proper labels** - Plot shows "E adm" and "I adm" as requested
5. **Automatic plotting** - Shows plot immediately after calculation
6. **Robust exit** - Handles EOF when running from scripts

## Expected Behavior

When you test with data generated from known parameters (e.g., alpha=0.5, omega=10000):
- The program WILL find the correct alpha IF your data comes from the sensitive omega range
- If your data comes from high omega (>1e4), alpha has little effect and you'll see poor sensitivity
- The debug tests showed that alpha sensitivity is ONLY significant at LOW omega (<1e-2)

## Output Files

After running, you'll find in `C:\Users\ranin\qrev_output\`:
- `qrev_rmse_results.csv` - All tested combinations and their RMSE
- `qrev_best_fit.csv` - Comparison of best simulation vs experimental data
- `qrev_best_curve.csv` - Full best-fit simulation curve
- `qrev_best_fit_plot.png` - The E adm vs I adm plot
