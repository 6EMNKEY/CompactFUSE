

if (!requireNamespace("ranger", quietly = TRUE)) {
  # Ranger package not installed, install it
  install.packages("ranger", repos = "https://cloud.r-project.org")
} else {
  # Ranger package is already installed
  print("ranger package is already installed.")
}


if (!requireNamespace("tidyverse", quietly = TRUE)) {
  # tidyverse package not installed, install it
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
} else {
  # tidyverse package is already installed
  print("tidyverse package is already installed.")
}

if (!requireNamespace("shinydashboard", quietly = TRUE)) {
  # shinydashboard package not installed, install it
  install.packages("shinydashboard", repos = "https://cloud.r-project.org")
} else {
  # shinydashboard package is already installed
  print("shinydashboard package is already installed.")
}

if (!requireNamespace("DT", quietly = TRUE)) {
  # DT package not installed, install it
  install.packages("DT", repos = "https://cloud.r-project.org")
} else {
  # DT package is already installed
  print("DT package is already installed.")
}