#======================================
# Loading the required packages 
#=======================================
# Load required libraries
library(sf)
library(haven)
library(grf)

#library(rgdal)
library(haven)

library(SpatialML)
library(GWmodel)      ## GW models
library(plyr)         ## Data management
library(sp)           ## Spatial Data management
library(spdep)        ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(raster)       ## spatial data
library(grid)         ## plot
library(gridExtra)    ## Multiple plot
library(ggplot2)      #  plotting
library(tidyverse)    # data 
library(SpatialML)    # Geographically weigted regression
library(dplyr)

library(grf)

library(tidyverse)
library(dplyr)
library(tmap)
#library(terra)

library(sf)
library(mapview)

library(rayshader) # 3d mapping
dhsgis<-st_read("fpc9.shp")
mapview(dhsgis)
length(unique(dhsgis$LATNUM))
dhsgis_complete<-filter(dhsgis, fvac>0)
dhs_df<-as.data.frame(dhsgis_complete)

#=======================================================
#SELECTED VARIABLES OF INTEREST AND STORED THEM IN A DF 
#=======================================================

suedhs_df <- dhs_df[, c("fvac","caseid","v005","v012","v021","v022","v023",
                        "v024","v025","v040","v106","v119","v123","v124",
                        "v125","v130","v131","v157","v158","v159",
                        "v169a","v171a","v190","v501","b2","b4",
                        "h1a", "LATNUM", "LONGNUM")]

##VARIABLE SELECTION##
#"bidx","v002","v003","v010","v034","v101","v102","v105a","
#v137","v150","v151","v153","v013","v120","v121","v155","v201","v394"","v504","bord","h1"

#vaccine variable removed from data set "hidx","h2","h3","h4","h5","h6","h7","h8","h9","h9a","h0","h10","h51","h52",
#"h53","h54","h55","h56","h57","h61","h62","h63","h64","h65","h66","s509y"

#======================================
# Mutating MEDIA and TRANSPORT variable. 
#======================================

suedhs_df <- suedhs_df %>%
  mutate(TRANSPORT = ifelse (
    v123== 1 | v124== 1 | v125== 1, 1,0
  ))



suedhs_df <- suedhs_df %>%
  mutate(MEDIA = case_when(
    v157 %in% c( 2, 3) ~ 1, # Received
    v158 %in% c( 2, 3) ~ 1,
    v159 %in% c( 2, 3) ~ 1,
    v169a %in% c( 1) ~ 1,
    v171a %in% c( 1) ~ 1,
    TRUE ~ 0        # Other cases (if any), assigning NA
  ))


#=====================================
# Divide data into train and test sets
#======================================
set.seed(123)  # for reproducibility
train_index <- sample(1:nrow(suedhs_df), size = 0.7 * nrow(suedhs_df))
train <- suedhs_df[train_index, ]
test <- suedhs_df[-train_index, ]

#======================
# DEFINE THE GWRF model
#======================

dhscoordsZ <- suedhs_df[, c("LATNUM", "LONGNUM")]  # XY coordinates

#=======================================
# Fit the Generalized Random Forest model               
#========================================
grf.model2 <- grf(fvac ~v005+v012+v021+v022+v023+
                    v024+v025+v040+v106+v119+TRANSPORT+v130+v131+MEDIA+
                    v190+v501+b2+b4+h1a,
                  dframe = suedhs_df,
                  coords = dhscoordsZ,
                  kernel = "adaptive",# Specify the kernel type
                  bw = 100,
                  weighted=TRUE  )

#=============
#model summary
#=============

summary(grf.model)
grf.model$LocalModelSummary
grf.model$Global.Model$importance

#======================
#ADD PREDICTIONS TO DF
#======================     # F= fvac preds

suedhs_df$f<-grf.model2[["Global.Model"]][["predictions"]]

suedhs_df$f

#==================================
#create point data of predictions
#==================================

fvacgis2<- suedhs_df[, c("LATNUM", "LONGNUM","f")]
head(fvacgis2)
fvacgis2<-st_as_sf(fvacgis2, coords = c("LONGNUM","LATNUM")) #convert to spatial
fvacgis2<-st_set_crs(fvacgis2, 4326) #define coordinate system

#write.csv(fvacgis,"fvacpred.csv") #write csv
#fvacgis<- read.csv("fvacpred.csv")
#fvacgis <- st_read("fvacpred.csv")

#================
#mapview fvacpred
#================
mapview(fvacgis, zcol="fvacpred")+
  #mapview (fvacgis,zcol="fvacpred")
  #mapview(counties, zcol="COUNTY", legend=FALSE)#preview
  
  # Plot
  mapview(fvacgis2, zcol = "f")
#==========================
# assign points to counties
#==========================
fvac_county<- st_join(fvacgis2, counties)

#==============================
#get median value for each county
#=============================
fvac_per_county2<- aggregate(fvac_county$f, list(fvac_county$COUNTY), FUN=median)
fvac_per_county2 <- fvac_per_county2 %>%
  dplyr::rename("COUNTY" = "Group.1")
fpc <- merge(counties, fvac_per_county2, by = "COUNTY")

mapview(fpc, zcol="x")

#=============================================
#Reclass FPC(MEDIANS)  to 1 or 2 BINARY OUTPUT
#=============================================
fpc3<- fpc %>%
  mutate(fvac_binary= case_when(
    x < 1.5 ~ "1", TRUE~"2"))
mapview(fpc3, zcol="fvac_binary")


#======================================================================
#MERGING DATA TO SHOW VARIABLE IMPORTANCE FOR EVERY VARIABLE PER COUNTY
#======================================================================

# merge 
fpc5<-fvac_county

fpc4<-(grf.model2[["Local.Variable.Importance"]])
length(unique(fpc4$v005))

fpc5$v005<-grf.model2[["Local.Variable.Importance"]][["v005"]]
fpc6<- merge(fpc5,fpc4, by="v005")

fpc7<- fpc6 %>% 
  group_by(COUNTY) %>%
  summarize_all(mean, na.rm=TRUE)

fpc8<- st_join(counties, fpc7)

fpc9<- subset(fpc8, select = -c(OBJECTID.x,  AREA.x,PERIMETER.x,COUNTY3_.x,v101.x,COUNTY.x,Shape_Leng.x,Shape_Area.x, Shape_Leng.y,Shape_Area.y,AREA.y, PERIMETER.y, COUNTY3_.y,OBJECTID.y,v101.y))
colnames(fpc9)
st_write(fpc9,"fpc9.shp")

#========================================================================================
#CREATING A SHINY DASHBOARD TO SHOW THE VAR IMPORTANCE OF EACH VARIABLE FOR EACH COUNTY
#=========================================================================================

library(shiny)
library(leaflet)
library(shinythemes)
library(sf)
library(rsconnect)


library(rsconnect)

## Adjust maximum bundle size limit
options(rsconnect.max.bundle.size = 2^31 - 1)

# Deploy app
rsconnect::deployApp()

# Load the shapefile or data
plotts <- st_read("fpc9.shp")


# Rename the column names of plotts to desired names
predshp2 <- c("COUNTY","Women Indv Samp Weight", "Predicted vaccination Medians","Respondents current age",
              "Pri Sampling Unit","Sample Strata for Errors","Stratification in sample design",
              "Region","Type of place of Residence","Cluster altitude in meters","Highest Education level",
              "Electricity","Transport","Religion","Ethnicity","Media","Wealth index combined",
              "Marital Status","Year of Birth","sex of child","Has Health card/ vaccination document","geometry")
names(plotts) <- predshp2

# Define UI
ui <- fluidPage(
  titlePanel(
    div(style = "font-weight: bold; font-family: 'Times New Roman, serif; color: #000000;",
        "GWRF MODEL LOCAL VARIABLE IMPORTANCE"
    ),
    windowTitle = "GWRF Model Importance"
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("attribute", "Select Attribute:", choices = colnames(subset(plotts, select = -c(COUNTY))), width = "100%"),
      tags$hr(),
      helpText("Select an attribute from the dropdown to visualize its importance on the map."),
      tags$hr()
    ),
    mainPanel(
      leafletOutput("map", width = "100%", height = "500px")
    )
  ),
  # 2. Change the background color of the app
  tags$style(HTML("body {background-color:#B0AC5F;}")),
  theme = shinytheme("darkly")
)

#Valid themes are: cerulean, cosmo, cyborg, darkly, flatly, journal,
#lumen, paper, readable, sandstone, simplex, slate, spacelab, superhero, united, yeti.

# Define server logic
server <- function(input, output, session) {
  
  output$map <- renderLeaflet({
    pal <- colorNumeric(
      palette = "Spectral",
      domain = plotts[[input$attribute]]
    )
    
    leaflet(plotts) %>%
      addProviderTiles("Esri.WorldImagery") %>%
      setView(lng = 37.5, lat = 0, zoom = 5) %>%
      addPolygons(
        fillColor = ~pal(plotts[[input$attribute]]),
        fillOpacity = 0.7,
        color = "#BDBDC3",
        weight = 1,
        stroke = TRUE,
        highlight = highlightOptions(
          weight = 2,
          color = "black",
          fillOpacity = 0.7,
          bringToFront = TRUE
        ),
        label = ~paste(COUNTY, ": ", plotts[[input$attribute]], sep = ""),
        labelOptions = labelOptions(
          style = list("font-weight" = "normal", padding = "3px 8px"),
          textsize = "15px",
          direction = "auto"
        )
      ) %>%
      addLegend(
        pal = pal,
        values = ~plotts[[input$attribute]],
        title = input$attribute,
        opacity = 0.7,
        position = "bottomright"
      )
  })
}

# Run the application
shinyApp(ui = ui, server = server)

