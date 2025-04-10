## R shiny app to simulate epidemics of a mosquito-borne disease
## and test factors determine the epidemic dynamics
## For P8477, by Wan Yang

library(shiny)
library(deSolve)
# function for a simple mosquito borne disease
SIRMosVec = function(time, state, parms) {
  with(as.list(c(state, parms)), {
    # infection in humans
    dSH = vH - SH * r * (THM * IM) - muH * SH
    dIH = SH * r * (THM * IM) - gamma * IH - muH * IH
    # dRH = gamma * IH - muH * RH
    
    # infection in mosquitoes
    dSM = vM - SM * r * (TMH * IH) - muM * SM
    dIM = SM * r * (TMH * IH) - muM * IM
    list(c(dSH, dIH, dSM, dIM))
  })
}

# parameters/initial conditions
NH=1e7; 
IH = 1; SH=NH-IH;
muH=1/50/365; # human life span: 50 yr
vH=NH*muH;
TMH = 0.2; # prob infection from human to mosquito;
THM = 0.1; # prob infection from mosquito to human;
gamma=1/7; # infectious period: 7 days
muM = 1/7; # 1 week life span
b=.5; # number of bite per day;
r = b / NH; # bite rate per human per mosquito

## NOTE: to find the maximum prevalence, time step needs to be small enough
## otherwise, may miss the max
tm.step=7; num_yr=500;
times=seq(0,365*num_yr,by=tm.step);



ui <- fluidPage(h2('Impact of Mosquito/Human ratio on Vector-borne Disease Dynamics'),
                
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'm2h.ratio',h6('Mosquito/Human ratio:'),value=6,min=1,max=10,step=.1),
                               h6('Number bites per day: 0.5'),
                               h6('Prob transmission from human to mosquito: 0.2'),
                               h6('Prob transmission from mosquito to human: 0.1'),
                               h6('Infectious period: 7 days'),
                               h6('Life span for human: 50 yr'),
                               h6('Life span for mosquito: 7 days')),
                              
                              mainPanel(
                                plotOutput(outputId = 'plots',width = '100%', height = "550px")
                                )
                              )
                
                )

server <- function(input, output){
  
  output$plots=renderPlot({
    
    NM=NH * input$m2h.ratio
    IM=1; SM=NM-IM;
    vM=NM*muM;
    state = c(SH = SH, IH = 1,SM = SM,  IM = 1)
    parameters = c(muH = muH, muM = muM,
                   vH = vH, vM = vM, 
                   THM = THM, TMH = TMH, 
                   gamma = gamma, r = b / NH)
    
    sim=ode(y=state,times=times,func=SIRMosVec,parms = parameters)
    
    R0=b^2*TMH*THM*NM/(muM*(gamma+muH)*NH)
    
    # plot results
    par(mfrow=c(2,1),mar=c(3,3,1,1), cex=1, mgp=c(1.5,.5,0))
    plot(times,sim[,'SH']/NH*100,type='l', xaxt='n',
            lwd=1.5,col='green',lty=1, main='Susceptibility in Human Population', cex.main=1.2,
            ylab='% Susceptible',xlab='Time (years)',ylim = c(0,100))
    axis(1,at=seq(0,tail(times,1),length.out = num_yr/10+1),labels = seq(0,round(length(times)/365*tm.step,0),by=10))
    mtext(paste0('Mosquito/Human ratio: ', round(NM/NH,3)),side=3,outer = F,line=-1.3,cex=1.2)
    mtext(paste0('R0: ', round(R0,2)),side=3,outer = F,line=-2.5,cex=1.2)
    
    plot(times,sim[,'IH']/NH*100,type='l', xaxt='n',
            lwd=1.5,col='red',lty=1,main='Prevalence in Human Population',cex.main=1.2,
            ylab='% Infectious',xlab='Time (years)',ylim = c(0,20))
    axis(1,at=seq(0,tail(times,1),length.out = num_yr/10+1),labels = seq(0,round(length(times)/365*tm.step,0),by=10))
    mtext(paste0('Mosquito/Human ratio: ', round(NM/NH,3)),side=3,outer = F,line=-1.3,cex=1.2)
    mtext(paste0('R0: ', round(R0,2)),side=3,outer = F,line=-2.5,cex=1.2)
  })
  
}

shinyApp(ui=ui, server = server)