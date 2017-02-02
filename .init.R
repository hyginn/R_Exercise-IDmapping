# .init.R
# Functions to initialize this Exercise session
# Boris Steipe
# ====================================================================

# Create a local copy of myScript.R if required, and not been done yet.
if (! file.exists("myScript.R" && file.exists(".tmp.R"))) {
    file.copy(".tmp.R", "myScript.R")
}

source(".utilities.R")

file.edit("R_Exercise-IDmapping.R")

# [End]
