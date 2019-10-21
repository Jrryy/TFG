NECESSARY PACKAGES FOR THE PDF REPORT GENERATION: 

texlive-latex-base
texlive-fonts-recommended
texlive-latex-extra

TO EXECUTE:

Open a terminal in the TFG_workflow directory, then execute:

$ R

Inside of the R interpreter:

> packages.install('shiny') # In case your computer does not have the shiny package installed
> library(shiny)
> runApp(appDir=".", launch.browser=FALSE, port=8000)

Then open localhost:8000 in any browser.

The mixOmics package might have unsatisfied C and C++ library dependencies. If its installation fails, installing the missing packages is enough to solve this problem.
