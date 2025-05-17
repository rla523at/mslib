@ echo off

for /R "_src" %%F in (*.cpp *.h) do (
  clang-format -i --style=file "%%F"
)

for /R "TEST_mslib" %%F in (*.cpp *.h) do (
  clang-format -i --style=file "%%F"
)

pause