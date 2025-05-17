@echo off

robocopy "../_src" "../__dist" "*.h" /njh /njs /nfl /ndl
robocopy "../_Output" "../__dist" "*.lib" /njh /njs /nfl /ndl

echo run!!!

