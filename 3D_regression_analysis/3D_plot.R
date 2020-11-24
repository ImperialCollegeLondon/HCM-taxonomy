########################### 3D Plot for UKBB data #################################
###################################################################################
###################################################################################
rm(list = ls(all = TRUE))  
#start with an empty environment

install.packages("plotly")
install.packages("data.table")
install.packages("dplyr")
# install all packages required

library(plotly)
library(data.table)
library(dplyr)

Data<-read.table("output/Sarc/HCM_beta_WT.txt")
colnames(Data)<-c("x","y","z","beta")
Data_pvalues<-read.table("output/Sarc/HCM_pvaluesTFCE_WT.txt")
colnames(Data_pvalues)<-c("x","y","z","pvalues")

# Some basic computations of average beta coefficient and the significance area
pos<-which(Data_pvalues$pvalues < 0.05)
sig <- (length(pos)/nrow(Data))*100
significance_area<-format(round(sig, 0), nsmall = 0)
average_beta <- format(round(mean(Data$beta[pos]), 3), nsmall = 2)

# Plot
ax <- list(
  title = "",
  zeroline = F,
  showline = F,
  showticklabels = F,
  showgrid = F)
p<-plot_ly(Data, x = ~x , y = ~y , z = ~z, 
           marker = list(color = ~beta,colorbar=list(title='beta coefficient - (a.u.)',len = 0.55, x = 1.00, y = 0.8), 
                         colorscale = c('#1972A4', '#FF7070'),
                         cauto = F,
                         cmin = -max(abs(Data$beta))*4,
                         cmax = max(abs(Data$beta))*4,showscale = F)) %>%
  add_markers() %>%
  layout(title = paste("\nComputational modelling of cardiac geometry in HCM subjects\n with average beta coefficient = ",average_beta,
                       "and significant area = ", significance_area, "%\n."),
         scene = list(xaxis = ax,yaxis = ax, zaxis = ax, dragmode="turntable",
                      camera = list(eye = list(x = cos(3.3)*2, y = sin(3.3)*2, z= 0.23))))

p<- add_trace(p = p,
              z = Data[pos,3],
              x = Data[pos,1],
              y = Data[pos,2],type = "scatter3d", mode = "markers",inherit = T,
              marker = list(color =  Data[pos,4],colorbar=list(title='beta coefficient - (a.u.)',len = 0.55, x = 1.13, y = 0.4),
                            colorscale = c('#1972A4', '#FF7070'),
                            cauto= F,
                            cmin = -max(abs(Data[pos,4])),
                            cmax = max(abs(Data[pos,4])), showscale = T))%>%
  layout(showlegend = F)
p
