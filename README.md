# shinyEVE
A shiny app for the Broad Institute group to vizualize EVE phylogenetic expression analysis


## installation notes
First git clone git@github.com:dylan-unlv/shinyEVE.git

Then unzip the data in the cloned directory, data.tar.gz

In order to get this working in its current form, make sure you have all the packages in src/necessary_packages.R installed into your R environment.

## running shinyEVE
Beginning in the same directory the repo is cloned to, use an R terminal by typing R on the command line

Then run shiny::runApp('src/shinyEVE') which should open a window with your default browser showing the shiny app

The app is still early in development, so give it a minute to start up. You'll know the app is loaded when the list of genes appears on the top row.

![plot](/images/test_image.png)
