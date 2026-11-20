@echo off
rem Double-click to run analysis 1 (control genes). Results are written into this folder.
cd /d "%~dp0"
where Rscript >nul 2>nul && Rscript 00_dump_genesets.R genesets_selected.gmt
where python >nul 2>nul && (python -I 01_control_genes.py) || (py -3 -I 01_control_genes.py)
echo.
echo Finished. You can close this window.
pause
