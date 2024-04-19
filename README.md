# shinyEVE
A shiny app for the Broad Institute to vizualize EVE phylogenetic expression analysis


## installation notes
First git clone git@github.com:dylan-unlv/shinyEVE.git

Then download the zipped data and store them in the same directory under the directory data

In order to get this working in its current form, make sure you have all the packages in src/shinyEVE/ui.R installed into your R environment.

You'll also have to change the local paths in src/shinyEVE/server.R in the first 30 or so lines

## running shinyEVE
Beginning in the same directory the repo is cloned to, use an R terminal by typing R on the command line

Then run shiny::runApp('src/shinyEVE') which should open a window with your default browser showing the shiny app
