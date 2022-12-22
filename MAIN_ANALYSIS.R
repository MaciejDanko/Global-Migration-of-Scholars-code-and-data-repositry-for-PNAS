###############################################################################
# Source code for performing analyses and plotting results and data in Figure 2 
# of the paper:       
#                   "Global Migration of Scholars: 
#         Trends and Patterns Revealed by Bibliometric Data" 
#
# CITATION: XXXXXXXXXXXXXXXXXX
# DOI: XXXXXXXXXXXX
#
# Author of the code: Maciej J. Dańko
# ORCID: 0000-0002-7924-9022
# Institution: Max Planck Institute for Demographic Research, Rostock, Germany
# WWW: https://github.com/MaciejDanko
# WWW2: https://www.researchgate.net/profile/Maciej-Danko
# Email: danko@demogr.mpg.de; maciej.danko@gmail.com
################################################################################

rm(list=ls()) 
gc()

require(mgcv)
require(parallel)
require(data.table)

################################################################################
## PREPARING DATA                                                             ##
################################################################################

# Load the data
load('AGGREGATED_DATA.RDA')
colnames(DAT)
length(unique(DAT$countrycode))

################################################################################
# Filter data by a GDP type and a number of countries of the highest average 
# number of scientists.
#
# Although Iceland and Luxembourg have a relatively large number of scientists, 
# they are not included in the list of countries for analysis because their 
# populations are very small (<500,000)
filter_data <- function(DAT, GDPType, NumberOfCountries, 
                        remove_small_countries = TRUE, minsize = 0.5){
  # Inits and checks
  GDPType <- GDPType[1] 
  ValidType <- colnames(DAT)[grep('GDP',colnames(DAT))]
  if (!GDPType%in%ValidType) stop(paste('Not valid GDPType, plese select one of',
                                        paste(paste('"', ValidType, '"', sep=''),
                                              collapse = ' , ')))
  NumberOfCountries <- NumberOfCountries[1]
  
  if (NumberOfCountries > 100) NumberOfCountries <- 100
  if (NumberOfCountries < 2) NumberOfCountries <- 2
  
  select_countries <- function(DAT, NumberOfCountries, cntr_names_var) {
    MeanScholars <- tapply(DAT$padded_population_of_researchers, 
                           DAT[,cntr_names_var], 
                           mean, na.rm=TRUE)
    OrderedCountries <- names(MeanScholars)[order(MeanScholars, 
                                                  decreasing = TRUE)]
    OrderedCountries[seq_len(NumberOfCountries)]
  }
  
  # Sort 
  DAT <- DAT[order(paste(DAT$countrycode, DAT$year)),]
  
  # Filter the GDP type and other other filtering and tests
  DAT$gdp <- DAT[,GDPType]
  UC_presel <- select_countries(DAT, NumberOfCountries, 'countryname')
  DAT <- DAT[which(!is.na(DAT$gdp)),]
  DAT$loggdp <- log10(DAT$gdp)
  DAT <- DAT[,which(!grepl('GDP_', colnames(DAT)))]
  DAT$GDPType <- GDPType
  
  # Remove small countries. For 100 countries these are Iceland, Luxembourg 
  if(remove_small_countries){
    smallCountries <- tapply(DAT$total_population_size, 
                             DAT$countrycode, mean, na.rm=TRUE)/1e6
    #minimum size in millions
    smallCountries <- names(smallCountries)[smallCountries < minsize] 
    TD <- unique(DAT$countryname[ which(paste(DAT$countrycode) %in% smallCountries )])
    DAT <- DAT[which(!DAT$countryname %in% TD),]
  }
  
  # Limit the number of countries to those that meet the criteria
  SelectedCountries <- select_countries(DAT, NumberOfCountries, 'countryname')
  NoGDP_coutries <- UC_presel[which(!UC_presel %in% SelectedCountries)]
  Replacement_countries <- SelectedCountries[
    which(!SelectedCountries %in% UC_presel)]
  if (length(NoGDP_coutries)) 
    message(paste0('The following countries have a ',  
                   'relatively large number of scholars but they had to be removed 
                   because off no GDP (', GDPType, ') data or being too small:\n\t',
                   paste(NoGDP_coutries, collapse = ', '), '.\n',
                   'These countries were replaced by the next in line countries',
                   'with the largest number of scholars:\n\t',
                   paste(Replacement_countries, collapse = ', '), '.'))
  
  DAT <- DAT[which(DAT$countryname %in% SelectedCountries),]
  
  DAT$countryname <- droplevels(DAT$countryname)
  DAT$countrycode <- droplevels(DAT$countrycode)
  DAT
}

# A dataset used for the main analysis
Main_Data <- filter_data(DAT, 
                       GDPType = 'GDP_PCAP_PPP_Const2017', 
                       NumberOfCountries = 100)

# Save list of 100 countries needed for Figure 1
tmp <- Main_Data[,c('countrycode', 'countryname')][
  which(!duplicated(Main_Data$countrycode)),]
colnames(tmp) <- c('ISO 3', 'NAME')
tmp <- tmp[order(tmp$NAME),]

MeanScholars <- tapply(Main_Data$padded_population_of_researchers, 
                       Main_Data$countryname, 
                       mean, na.rm=TRUE)

MeanScholars <- MeanScholars[order(names(MeanScholars))]

if (any(paste(tmp$NAME) != names(MeanScholars))) stop()

tmp$`AVG SCHOLARS` <- as.numeric(MeanScholars)
tmp <- tmp[order(tmp$`AVG SCHOLARS`, decreasing = TRUE),]

write.csv(tmp, paste('COUNTRIES_LIST_', length(unique(Main_Data$countrycode)),
                     '.CSV', sep = ''), row.names = FALSE)

if (FALSE) { 
  # Alternative GDP measures used for the sensitivity analysis
  Alt_Data_1 <- filter_data(DAT, 
                          GDPType = 'GDP_PCAP_PPP_Current', 
                          NumberOfCountries = 100)
  length(unique(Alt_Data_1$countrycode))
  
  Alt_Data_2 <- filter_data(DAT, 
                          GDPType = 'GDP_PCAP_Const2010', 
                          NumberOfCountries = 100)
  length(unique(Alt_Data_2$countrycode))
}


################################################################################
## FITTING MODEL(S)                                                           ##
################################################################################

# Main model formula used in the paper
Main_Formula <- number_of_outmigrations ~
  s(countrycode, bs = "re") + 
  s(year, bs = 'ps') + 
  s(loggdp, bs = 'ps') + 
  s(countrycode, loggdp, bs = "re") +
  offset(logexposures)

# Simple model formula used in the paper
Simple_Formula <- number_of_outmigrations ~
  s(year, bs = 'ps') + 
  s(loggdp, bs = 'ps') + 
  offset(logexposures)

# Main model used in the paper
Main_Model <- mgcv::gam(Main_Formula,
                        data = Main_Data,
                        family = nb(),
                        # remove "control" on error 
                        control = gam.control(nthreads = 
                                                parallel::detectCores() - 1, 
                                              epsilon = 1e-08, 
                                              maxit = 500,
                                              efs.tol = 0.025, 
                                              edge.correct = TRUE),
                        select = FALSE,
                        method = 'REML')

# Simple model used in the paper
Simple_Model <- mgcv::gam(Simple_Formula,
                          data = Main_Data,
                          family = nb(),
                          select = FALSE,
                          method = 'REML')

# Compare models
AIC(Main_Model, Simple_Model)

# Smooth term of year
mgcv:::plot.gam(Main_Model, pages = 0, select=2, ylab='s(year)', ylim=c(-0.5, 0.5))

# Smooth term of log GDP
mgcv:::plot.gam(Main_Model, pages = 0, select=3, ylab='s(log GDP)', xlab='log GDP')


################################################################################
## PLOT RESULTS                                                               ##
################################################################################

# Fast merging data frames function
fast.merge.df<-function(DF1, DF2, by, all.x = TRUE, all.y = TRUE){
  DT1 <- data.table :: data.table(DF1, key = by, stringsAsFactors = FALSE)
  DT2 <- data.table :: data.table(DF2, key = by, stringsAsFactors = FALSE)
  data.frame(data.table ::: merge.data.table(DT1, DT2, 
                                           all.x = all.x, all.y = all.y),
             stringsAsFactors = FALSE)
}

# Main plotting function
plot_data_and_fit <- function(Main_Model, 
                              Simple_Model,
                              Data,
                              show_country_specific_data = TRUE, 
                              show_country_specific_fits = TRUE, 
                              plot_simple_model_fit = TRUE,
                              ylim = c(5,200), 
                              logGDP_class_size = 0.25,
                              ax_font = 1,
                              xnam = 'GDP per capita in 1,000s of US dollars', 
                              ynam = 'Emigration rate per 1,000 scholars'){
  
  rowMeans_ <- function(x, ...){
    x <- as.matrix(x)
    if(any(dim(x) == 1)) x else rowMeans(x, ...)
  }
  
  rowSums_<-function(x, ...){
    x <- as.matrix(x)
    if(any(dim(x) == 1)) x else rowSums(x, ...)
  }
  
  inf2NA <- function(x) {x[is.infinite(x)] <- NA; x}
  
  ##############################################################################
  # Initial calculations and settings
  ##############################################################################
  
  YEAR <- unique(Data$year)
  CNTR <- levels(droplevels(Data$countrycode))
  CNTRPAL <- adjustcolor(rainbow(length(CNTR) + 1)[-1], alpha.f =  0.3)
  DATN <- length(unique(Data$countrycode))
  OutMigRate <- Main_Data$number_of_outmigrations/
    Main_Data$padded_population_of_researchers
  
  ylim <- log(ylim / 1000)
  minX <- floor(min(Main_Data$loggdp) / logGDP_class_size) * logGDP_class_size
  maxX <- ceiling(max(Main_Data$loggdp) / logGDP_class_size) * logGDP_class_size
  LogGDPBreaks <- seq(minX, maxX, logGDP_class_size)
  ClassLogGDP <- cut(Main_Data$loggdp, breaks = LogGDPBreaks)
  
  x2pos <- function(x) (x - LogGDPBreaks[1]) * (1 / logGDP_class_size) + 0.5
  PredX <- seq(LogGDPBreaks[1], LogGDPBreaks[length(LogGDPBreaks)], 0.01)
  xx <- x2pos(PredX)
  
  if (length(ylim) == 0){
    ylim <- range(inf2NA(log(yrat)), na.rm = TRUE)
    ylim[2] <- max(ylim[2], -1.3)
    ylim[1] <- min(ylim[1], -7)
  }
  
  # Preparing a data frame for country-specific model predictions
  if (show_country_specific_fits){
    WWW <- expand.grid( CNTR = CNTR, YEAR = YEAR)
    WWW$ind <- paste(WWW$CNTR, WWW$YEAR, sep = '_')
    WWW <- fast.merge.df(
      data.frame(Exposures = Main_Data$padded_population_of_researchers,
                 ind = paste(Main_Data$countrycode, Main_Data$year, sep = '_')),
      WWW, by = 'ind')
    WWW <- fast.merge.df(
      data.frame(LogGDP = Main_Data$loggdp,
                 ind = paste(Main_Data$countrycode, Main_Data$year, sep = '_')),
      WWW, by = 'ind')
  }
  
  ##############################################################################
  # Calculation of year-averaged predicted log emigration rates described in 
  # supplementary materials (eq. 1)
  ##############################################################################
  
  # Main model
  sterms<-rownames(summary(Main_Model)$s.tab)
  excludeterms<-sterms[grep('country', sterms)]
  excluZZ<-sapply(YEAR,
                  function(k) 
                    predict.gam(Main_Model, 
                                data.frame(loggdp = PredX, 
                                           year = k, 
                                    # can be any country, it is excluded anyway
                                           countrycode = 'POL', 
                                           logexposures = 0), 
                                exclude = excludeterms,
                                type = 'link',
                                se.fit = FALSE))
  
  year_averaged_predicted_log_emigration_rates <- rowMeans_(excluZZ)
  
  # Simple model
  ZZ0<-sapply(YEAR,function(k) 
    predict(Simple_Model,
            data.frame(loggdp = PredX, 
                       year = k, 
                       logexposures = 0),
            type = 'link',
            se.fit = FALSE))
  
  year_averaged_predicted_log_emigration_rates_simple <- rowMeans_(ZZ0)
  
  # Main model, country-specific fits
  if (show_country_specific_fits){
    ZZ <- lapply(1 : nrow(WWW), function(k) 
      predict.gam(Main_Model, 
                  data.frame(loggdp = PredX, 
                             year = WWW$YEAR[k], 
                             countrycode = WWW$CNTR[k], 
                             logexposures = 0),
                  type = 'link',
                  se.fit = FALSE))
    ZZ<-sapply(ZZ,function(k) k)
    CNTRV <- sapply(CNTR, function(k) 
      rowMeans_(ZZ[,which(paste(WWW$CNTR) == k)]))
    CNTRR <- sapply(CNTR, function(k) 
      range(Main_Data$loggdp[Main_Data$countrycode == k]))
  }
  
  ##############################################################################
  # plotting figure 2
  ##############################################################################
  
  # plot box-plots, y-axis, labels, title, etc
  plot(ClassLogGDP, inf2NA(log(OutMigRate)), ylim = ylim, ylab = '', xlab = '',
       col = adjustcolor('gray', 0.5), border = adjustcolor('black', 0.5),
       las = 3, axes = FALSE)
  ay <- c(0.00001, 0.000025, 0.00005, 0.0001, 0.00025, 0.0005, 0.001, 0.0025,
        0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1)
  axis(2, log(ay), ay * 1000, las = 1, font = ax_font)
  box();
  
  mtext(xnam, 1, 2.4, font = ax_font)
  mtext(ynam, 2, 2.4, font = ax_font)
  if (length(YEAR) > 1) mtext(paste(paste(min(YEAR), '-', max(YEAR)), ', ', DATN,
                        ' Countries', sep = ''), 3, 0.1, font = ax_font) else 
                        mtext(paste(YEAR, ', ', DATN, ' Countries',
                        sep = ''), 3, 0.1, font = ax_font)
  
  # plot country-specific fits
  if (show_country_specific_fits){
    for (k in 1 : length(CNTR)) {
      plo <- approx(xx, CNTRV[,k], x2pos(unique(seq(CNTRR[1,k],
                                                CNTRR[2,k],
                                                length.out=100))))
      if (length(plo$x) > 1) typee <- 'l' else typee <- 'p'
      lines(plo, col = CNTRPAL[k], type = typee, lwd = 2)
    }
  }
  
  # plot country-specific data points
  if (show_country_specific_data)
    lines(x2pos(Main_Data$loggdp), log(OutMigRate), type = 'p',
          col = CNTRPAL[as.numeric(Main_Data$countrycode)], pch = 20, cex = 0.5)
  
  # plot main model marginal fit
  lines(xx,year_averaged_predicted_log_emigration_rates, 
        col = 1, lwd = 3.5, type = 'l')
  lines(xx,year_averaged_predicted_log_emigration_rates,
        col = adjustcolor('darkorange', alpha.f = 0.8),
        lwd = 3, type = 'l')
  
  # plot simple model marginal fit
  if(plot_simple_model_fit) {
    lines(xx,year_averaged_predicted_log_emigration_rates_simple,
          col = 1, lwd = 3.5, type = 'l')
    lines(xx,year_averaged_predicted_log_emigration_rates_simple,
          col = adjustcolor(4, alpha.f = 0.8), lwd = 3, type = 'l')
  }
  
  # plot x-axis
  gdplab <- c(50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 
              100000, 250000)
  xxl <- x2pos(log10(gdplab)) 
  axis(1, xxl, labels = gdplab / 1000, font = ax_font)
  abline(v = x2pos(seq(1,6,0.25)), col = 'gray')
  
  # plot legend if simple model is also plotted
  if(plot_simple_model_fit) {
    legend('bottomleft', c('Without country random-effects',
                        'With country random-effects'),
           lwd = 3.5, bty = 'n')
    legend('bottomleft', c('Without country random-effects',
                        'With country random-effects'),
           col=c(adjustcolor(4, alpha.f = 0.8),
                 adjustcolor('darkorange', alpha.f = 0.8)),
           lwd = 3, bty = 'n')
  }
}

# Save the plot
pdf('./FIGURES/Figure_2.pdf', width=9.73, height=7)
par(mfrow = c(1, 1))
par(mar = c(3.4, 3.4, 1, 0.5), oma=c(0,0,0,0))
plot_data_and_fit(Main_Model, 
                  Simple_Model,
                  Main_Data,
                  show_country_specific_data = TRUE, 
                  show_country_specific_fits = TRUE, 
                  plot_simple_model_fit = TRUE,
                  ylim = c(4.9916, 165.2989), 
                  logGDP_class_size = 0.25,
                  ax_font = 2,
                  xnam = 'GDP per capita in 1,000s of US dollars', 
                  ynam = 'Emigration rate per 1,000 scholars')
dev.off()
