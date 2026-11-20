@echo off
rem Double-click to run step 2 of the generalization analysis. Results go to the "results" subfolder.
rem Expected run time: about 1-2 hours. For a 5-minute test run first, use RUN_generalize_quick.bat
cd /d "%~dp0"
if not exist genesets_all.gmt (where Rscript >nul 2>nul && Rscript 00_dump_genesets_all.R genesets_all.gmt)
where python >nul 2>nul && (python -I 01_generalize.py) || (py -3 -I 01_generalize.py)
echo.
echo Finished. You can close this window.
pause
