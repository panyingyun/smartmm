set source_mesh=%1
set target_mesh=%2
OT.exe -nearest -source %source_mesh%.m -target %target_mesh%.m -step_length 0.05
