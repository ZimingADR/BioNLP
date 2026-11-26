@echo off
title LaTeX Auto-Builder
chcp 65001 >nul
cd /d "%~dp0\paper_bibm\BIBM-Latex"
python auto_build.py
pause
