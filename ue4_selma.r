#Package ? utiliser 
library(packHV)
library(epiDisplay)
library(survival)
library(prettyR)

library(dplyr)
library(DT)
library(shiny)

library(ggplot2)
library(plotly)

library(mlbench)
library(rattle)

library(rpart)
library(rpart.plot)
library(RColorBrewer)

library(randomForest)

#Importer les donn?es patients### 
#dataLung <-file.choose()
dataLung = "Lung3.metadata.csv"


#lecture de la table dataLung
datap <- read.csv(file = paste0(dataLung), sep=";", header=T)
datap$characteristics.tag.tumor.size.maximumdiameter <- as.numeric(as.character(datap$characteristics.tag.tumor.size.maximumdiameter))

geo <- read.csv(file = 'data.csv', sep=";", header=T)
top10_matrix <- read.csv(file = 'top10_matrix.csv', sep=",", header=T)

#Choisir les variables les plus pertinentes 
dttp = datap[,c("title","source.location",
                "characteristics.tag.gender",
                "characteristics.tag.histology",
                "characteristics.tag.tumor.size.maximumdiameter",
                "characteristics.tag.stage.primary.tumor",
                "characteristics.tag.stage.nodes",
                "characteristics.tag.stage.mets",
                "characteristics.tag.grade")]

# Renommer les variables 
names(dttp) <- c("patient","Localisation","Gender","Histology",
                 "tumor.size","stage.T","stage.N",
                 "stage.M","grade")




# On regroupe en 2 groupes histologiques
dttp$Histology <- ifelse(
    dttp$Histology == "Adenocarcinoma With Mucinous Features"
    | dttp$Histology == "Adenocarcinoma With Papillary Features"
    | dttp$Histology == "Adenocarcinoma, Acinar Type"
    | dttp$Histology == "Adenocarcinoma, Bronchiolo-alveolar Features"
    | dttp$Histology == "Adenocarcinoma, Mixed"
    | dttp$Histology == "Adenocarcinoma, Mucinous With Bronchiolo-alveolar Features"
    | dttp$Histology == "Adenocarcinoma, NOS"
    | dttp$Histology == "Adenocarcinoma, Papillary, NOS"
    | dttp$Histology == "Papillary Type AND Adenocarcinoma, Bronchiolo-alveolar Features"
    
    , "Adenocarcinoma"
    ,  ifelse(dttp$Histology == "Squamous Cell Carcinoma, NOS", "Squamous_Carcinoma_NOS", "Autres")
  )



graphique <- function(donnees){
  
  tmpTab <- table(subset(dttp, select=c(donnees)))
  
  df <- data.frame(
    group = names(tmpTab),
    value = c(tmpTab)
  )
  
  ggplot(df, aes(x=group, y=value))+
    geom_bar(stat = "identity", aes(fill = group))
  

}


#Analyse descriptive

# 1-Localisation: Variable qualitative
summary(dttp$Localisation)
pie(table(dttp$Localisation))

# 2-Gender : variable qualitative
summary(dttp$Gender)
plot(dttp$Gender,main="Gender",col="pink")

# 3- Histology: Variable qualitative
#Regrouper les moins exprim?s en categorie "autres"
table(datap$characteristics.tag.histology)
summary(dttp$Histology)
#levels(dttp$Histology)<-c("autres", "autres", "autres", 
#                                          "Adenocarcinoma, Bronchiolo-alveolar Features", 
#                                          "autres","autres", "Adenocarcinoma, NOS",
#                                          "Adenocarcinoma, Papillary, NOS", "autres",
#                                          "Non-Small Cell", "autres","autres","autres",
#                                          "Squamous Cell Carcinoma, NOS","autres", "autres",
#                                          "autres")
#levels(dttp$Histology)
#describe(dttp$Histology)
#plot(dttp$Histology,main="characteristics.tag.histology",col="green")

# 4Tumor size : variable quantitative

summary(dttp$tumor.size)
table(dttp$tumor.size)
describe(table(dttp$tumor.size))
barplot(table(dttp$tumor.size), 
        col = "purple", 
        border = "white",
        main ="Tumor.size",
        xlab = "size")

# 5- stage.T: variable 
# Recodage de la variable stage.T pour grouper les sous categories
#Stage de la tumeur
levels(dttp$stage.T)<-c("pT1","pT1","pT1","pT2","pT2","pT2",
                                    "pT2","pT3","pT4","pTx", "pTx")
levels(dttp$stage.T)
table(dttp$stage.T)
pie(table(dttp$stage.T))

# 6-stage.N: variable qualitative
summary(dttp$stage.N)
plot(dttp$stage.N,main="characteristics.tag.stage.N",col="blue")

# 7-Stage.mets : variable qualitative
summary(dttp$stage.M)
plot(dttp$stage.M,main="characteristics.tag.stage.mets",col="magenta")

#8- characteristics.tag.grade: variable qualitative 
summary(dttp$grade)
plot(dttp$grade,main="characteristics.tag.grade",col="red")



#GSE58661_matrix = read.csv('GSE58661_series_matrix.txt', sep="\t", header=T, comment.char="!")
#top_10 =  geo$ID[1:10]
#top10_GSE58661_matrix = GSE58661_matrix[GSE58661_matrix$ID_REF %in% top_10, ]
#write.csv(top10_GSE58661_matrix, file = "top10_matrix.csv")


tree <- function( choix ){
  #matrice des 10 expressions significatives
  tab_top = top10_matrix
  
  #transpose
  tab_top_transpose = data.frame(t(tab_top))
  
  #renomme colonnes
  mod = as.character(t(tab_top[2]))
  mod = gsub( "-", "_",mod)
  names(tab_top_transpose) = mod
  
  tab_top_transpose = tab_top_transpose[-1,]
  tab_top_transpose = tab_top_transpose[-1,]
  tab_top_transpose$Histology = dttp$Histology
  
  DataNew = tab_top_transpose
  DataNew$Histology = as.factor(DataNew$Histology)
  
  
  # Construction des echantillons d'apprentissage et des echantillons de test
  set.seed(111) # initialisation du generateur
  # Extraction des ?chantillons
  test.ratio=.2 # part de l'?chantillon test
  npop=nrow(DataNew) # nombre de lignes dans les donn?es
  nvar=ncol(DataNew) # nombre de colonnes
  # taille de l'?chantillon test
  ntest=ceiling(npop*test.ratio)
  # indices de l'?chantillon test
  testi=sample(1:npop,ntest)
  # indices de l'?chantillon d'apprentissage
  appri=setdiff(1:npop,testi)
  # construction de l'?chantillon d'apprentissage
  datapq=DataNew[appri,]
  # construction de l'?chantillon test
  datestq=DataNew[testi,]
  # summary(datapq)
  
  for (i in 1:(10)){datapq[,i]=as.numeric(as.character(datapq[,i]))}
  for (i in 1:(10)){datestq[,i]=as.numeric(as.character(datestq[,i]))}
  
  # Exemple de mod?lisation avec un arbre de d?cision
  if (choix == "TREE"){
    fitq.tree=rpart(Histology~.,data=datapq,
                    parms=list(split= "application" ),method="class",
                    control=rpart.control(minsplit= 1 ,cp=0.05))
    fancyRpartPlot(fitq.tree)
  }
  
  if (choix == "forest"){
    # For?t al?atoire
    # Avec les seules variables quantitatives :
    fit=randomForest(Histology~.,data=datapq,do.trace=50,
                     importance=TRUE,norm.vote=FALSE)
    # Importance de chaque variable
    print(round(fit$importance, 2))
    print(table(predict(fit,datestq),datestq$Histology))
  }
}



### shiny
ui_0 <- fluidPage(

  sidebarLayout(
    sidebarPanel (
      
      selectInput(
        "patient", 
        "Choix patient : ",
        dttp$patient
      ),
      
      selectInput(
        "donnees", 
        "Graphique : ",
        colnames(dttp[-1])
        
      )
      
    ),
    
    mainPanel (
      tabsetPanel(
        tabPanel("Patients",DT::dataTableOutput('patients_out')),
        tabPanel("Patient",DT::dataTableOutput('patient_out')),
        
        tabPanel("Graphique", plotlyOutput('graphiques')),
        
        tabPanel("Geo", DT::dataTableOutput('geo_out')),
        
        tabPanel("TREE", plotOutput('arbre_out')),
        tabPanel("Random Forest", verbatimTextOutput('foret_out'))
        
      )
      
    )
    
  )
)  


server_0 <- function(input, output) {
  
  output$patients_out  <- DT::renderDataTable(
    dttp
  )
  
  output$patient_out <- DT::renderDataTable(
    dttp %>% 
      filter(patient == input$patient)
  )

  #Graphique :
  output$graphiques <- renderPlotly(
    graphique(as.character(input$donnees))
  )
  
  #geo2r :
  output$geo_out <- DT::renderDataTable(
    geo
  )
  
  #tree
  output$arbre_out <- renderPlot(
    tree("TREE")
  )
  
  output$foret_out <- renderPrint(
    tree("forest")
  )
  
}

shinyApp(ui = ui_0, server = server_0)
