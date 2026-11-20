@echo off
rem Double-click to run a QUICK TEST of step 2. Results go to the "results" subfolder.
rem Expected run time: about 1-2 hours. For a 5-minute test run first, use RUN_generalize_quick.bat
cd /d "%~dp0"
if not exist genesets_all.gmt (where Rscript >nul 2>nul && Rscript 00_dump_genesets_all.R genesets_all.gmt)
where python >nul 2>nul && (python -I 01_generalize.py --quick) || (py -3 -I 01_generalize.py --quick)
echo.
echo Finished. You can close this window.
pause
