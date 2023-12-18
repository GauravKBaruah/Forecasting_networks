library(statmod)

cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)



network_structure<-function(Na,Np, g ){
  new_g<- matrix(0, nrow= length(Np), ncol =length(Na) )
  
  Na[which(Na < 5e-1)]<-0
  Np[which(Np < 5e-1)]<-0
  
  for(i in 1:length(Np)){
    for(j in 1:length(Na)){
      new_g[i,j]<-g[i,j]*Na[j]*Np[i]    
      
    }
  }
  new_g[which(new_g > 0)]<-1
  
  return(new_g) 
}


#conversion to a matrix
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(unname(d))
  dat[dat > 0] = 1
  dat[dat < 0] = 0
  dat<-apply(dat,2,as.numeric)
  return(dat)}



#function for estimating trait matching for a network
trait.matching<-function(mA,mP,adj.mat_1,gamma){
  tm<-numeric()
  for(i in 1:nrow(adj.mat_1)){
    tm[i] <- mean(adj.mat_1[i,]*exp(-(mA-mP[i])^2)/gamma)
    
  }
  return(tm=mean(tm))
}



#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.animals<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.animal){
  
  
  temp2<-dat2<-x2<-x3<-array(dim=c(points))
  if(mat == 0){
    return(list(G= 0, B = 0))
  }
  else if(mat == 1){
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sigma$sa)$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      
      
      f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
      
      temp2[i]<- sum(np*(mut.strength/degree.animal)*f/(1+h*np*(mut.strength/degree.animal)*f)*w2*w1[i])
      
      
      
      dat2[i]<- sum( ((z1[i]-m$ma)*f*(mut.strength/degree.animal)*np/(1+h*np*(mut.strength/degree.animal)*f) )*w2*w1[i])
      
      
      
    }
    G = sum(temp2)
    B = sum(dat2) 
    
    
    
    
    return(list(G= G, B = B))
  }
}

#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.plants<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant){
  
  temp2<-dat2<-x3<-x4<-array(dim=c(points))
  
  if(mat == 0){
    
    return(list(G= 0, 
                B = 0))
    
  }
  else if (mat==1){
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sigma$sp)$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      
      f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
      
      temp2[i]<- sum(na*(mut.strength/degree.plant)*f/(1+h*na*(mut.strength/degree.plant)*f)*w2*w1[i])
      
      
      
      dat2[i]<- sum( ((z1[i]-m$mp)*f*(mut.strength/degree.plant)*na/(1+h*na*(mut.strength/degree.plant)*f) )*w2*w1[i])
      
      # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
    }
    
    G= sum(temp2)
    B = sum(dat2)
    
  } 
  
  
  
  return(list(G= G, 
              B = B))
  
  
}

cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)

eqs_euler <- function(time, y, pars) {
  A <- dim(pars$matrix)[2]  
  P <-dim(pars$matrix)[1]
  Np<- muP<-matrix(0, nrow=time, ncol=P)
  Na<- muA<-matrix(0, nrow =time, ncol=A)
  Na[1,]<-y[1:A]
  Np[1,]<-y[(A+1):(A+P)]
  muA[1,]<-y[(A+P+1):(A+P+A)]
  muP[1,]<-y[(A+P+A+1):(A+P+A+P)]
  ## define g, where g[i] is the selection pressure on species i from growth
  alpha.a<-pars$Amatrix # alpha matrix
  alpha.p<-pars$Pmatrix
  s<-pars$sig
  #w<-pars$w
  aij<-bij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-matrix(0, nrow=P,ncol=A) 
  aj<-bj<-ai<-bi<-numeric()
  degree.animals<-degree.plants<-numeric()
  
  dt<-pars$dt
  
  #degree of plants and animals
  for(i in 1:P){
    degree.plants[i]<-sum(pars$matrix[i,])} # degree of plants
  for(j in 1:A){
    degree.animals[j]<-sum(pars$matrix[,j]) # degree of animals
  }
  
  temper<-Temp<-numeric()
  temper[1]<-pars$initial_temperature
  Temp<-temper[1]
  aj<-bj<-ai<-bi<-vi<-vj<-numeric()
  for (t in 1:(time-1)){
    
    temper[t+1] <- temper[t] + pars$rate * pars$dt #*temper[t]*pars$dt#+ rnorm(1,0,2)
    Temp[t+1] <- rnorm(1,temper[t+1],sd=pars$env_var)
    #growth rate of animals and plants as a function of current temperature locally
    ba<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+s[1:A]))*exp(-(Temp[t]- muA[t,])^2/(2*(pars$bw)^2+s[1:A])) - 0.2*exp(10^4*(1/(muA[t,]+273)- 1/(Temp[t]+273)))
    bp<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+s[(A+1):(A+P)]))*exp(-(Temp[t]- muP[t,])^2/(2*(pars$bw)^2+s[(A+1):(A+P)])) #- 0.2*exp(10^4*(1/(muP[t,]+273)- 1/(Temp[t]+273)))
    
    bar_ba<- pars$gi/(pars$bw)*exp(-(Temp[t]- muA[t,])^2/(2*(pars$bw)^2+s[1:A]))*s[1:A]*(pars$bw)*(Temp[t] -muA[t,])/((pars$bw)^2 + s[1:A])^1.5
    bar_bp<-pars$gi/(pars$bw)*exp(-(Temp[t]- muP[t,])^2/(2*(pars$bw)^2+s[(A+1):(A+P)]))*s[(A+1):(A+P)]*(pars$bw)*(Temp[t] - muP[t,])/((pars$bw)^2 + s[(A+1):(A+P)])^1.5 
    
    for(r in 1:A){
      for(l in 1:P){
        #
        m.temp<-list(ma=muA[t,r],mp=muP[t,l])
        sigma1<-list(sa=s[r],sp=s[(A)+l])
        
        
        temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np[t,l],na=Na[t,r],
                                mut.strength=pars$mut.strength, points=5,
                                mat=pars$matrix[l,r],
                                degree.animal = degree.animals[r])
        aij[r,l]<-temp1$G
        bij[r,l]<-temp1$B
        
      }
      ai[r]<-sum(aij[r,])
      bi[r]<-sum(bij[r,])
      
    }
    for(k in 1:P){
      for(m in 1:A){
        m2.temp<-list(ma=muA[t,m],mp=muP[t,k])
        sigma2<-list(sa=(s[m]),sp=(s[(A+k)]))
        temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np[t,k],na=Na[t,m],
                               mut.strength=pars$mut.strength,
                               points=5,mat=pars$matrix[k,m], 
                               degree.plant =degree.plants[k])
        aji[k,m] <-temp2$G
        bji[k,m]<-temp2$B
        
      }
      aj[k]<-sum(aji[k,])
      bj[k]<-sum(bji[k,])
      
    }
    
    
    
    Na[t+1, ]<- ifelse(Na[t, ] < 1e-6, 0, Na[t,] + Na[t,]*(ba-alpha.a%*%Na[t,]+ai)*pars$dt + rnorm(A,0,0.1))# population dynamics
    Np[t+1, ]<- ifelse(Np[t, ] < 1e-6, 0, Np[t,] + Np[t,]*(bp-alpha.p%*%Np[t,]+aj)*pars$dt + rnorm(P,0,0.1))# population dynamics
    muA[t+1, ]<- muA[t,]+ pars$h2[1:A]*(bar_ba+ bij%*%Np[t,])*pars$dt                                         # mean trait dynamics
    muP[t+1, ]<-  muP[t,]+ pars$h2[1:P]*(bar_bp+ bji%*%Na[t,])*pars$dt                                        # mean trait dynamics
    
    Na[t+1,which(Na[t+1,] < 1e-4)]<-0
    Np[t+1,which(Np[t+1,] < 1e-4)]<-0

  }
  
  return(list(Na=Na, Np=Np,muA=muA,muP=muP,Temp=Temp))
}


## Organize simulation results into tidy table (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma
organize_results <- function(sol, pars) {
  S <- length(pars$sigma) ## number of species
  A<-dim(pars$matrix)[2] # no. of animals
  P<-dim(pars$matrix)[1] # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  # temp<- temp %>% filter(time >= pars$cutoff.time)
  names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)
  names(temp)[1] <- "time"
  names(temp)[(A+1+1):(A+P+1)] <- paste0("P_",  A+1:P) ## name trait mean columns
  names(temp)[(A+P+2):(A+P+A+1)] <- paste0("ma_",1:A)
  names(temp)[(A+P+A+2):(A+P+A+P+1)]<-paste0("mp_",  A+1:P)
  temp <- temp %>%
    gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
    separate(variable, c("type", "species"), sep="_") %>%
    # spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species, v) %>% ## rearrange columns
    mutate(species=as.integer(species), Temp=pars$Temp) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}


#params : list of parameter values
#ic = initial conditions
#tmax= time for simulations

cluster_run<-function(params, ic, tmax ){
  
  
  
  sol<-ode(func=eqs_noise, y=ic, parms=params, times=seq(0, tmax, by=tmax/500)) %>% ## solve ODEs
    organize_results(params) 
  
  A<-params$A  # no. of animal species
  P<- params$P # no. of plant species
  sol<-  sol %>% filter(time == tmax)
  
  Na<-sol %>% filter(type =="N")
  Np<-sol %>% filter(type =="P")
  ma<-sol %>% filter(type == "ma")
  mp<-sol %>% filter(type == "mp")
  
  
  richness_animal<- length( which(Na$v>1e-3 ) ==TRUE)
  richness_plant<- length( which(Np$v>1e-3 ) ==TRUE)
  
  trait_match<-matrix(0,nrow=A,ncol=P)
  
  if(params$var== "high"){ s <- 0.09}else{s<-0.009}
  
  muA <- (ma %>% filter(time == tmax))$v
  muP <-(mp %>% filter(time == tmax))$v
  
  for(i in 1: A){
    
    trait_match[i,] <- exp(-(muA[i]- muP)^2/(2*s+1)) 
  }
  
  trait_matching<-mean(trait_match,na.rm=T)
  
  
  
  Na$v[Na$v<0] <- 0
  Np$v[Np$v<0] <- 0
  pa<- Na$v/sum(Na$v)
  pp<-Np$v/sum(Np$v)
  #sum(pa^2)^(1/(1-2))
  richness<- richness_plant+richness_animal
  diversity<- (sum(pa^2)^(1/(1-2))+ sum(pp^2)^(1/(1-2)))/2
  biomass<- sum(Na$v)+sum(Np$v)
  
  
  
  #  
  output_network_data<- data.frame(richness=richness,
                                   diversity=diversity,
                                   biomass=biomass,
                                   trait_matching=trait_matching,
                                   individual_variation=params$var,
                                   h2= params$h2,
                                   Temperature=params$Temp,
                                   web=params$web,
                                   mut.strength=params$mut.strength,
                                   nestedness=params$nestedness,
                                   connectance=params$C,
                                   network_size= params$network_size)
  
  
  
  
  species_level_data <-as.data.frame(cbind(
    rep(seq(1,(nrow(params$matrix)+ncol(params$matrix)),1)), #number of species
    as.numeric(c( (Na %>% filter(time == tmax))$v, (Np %>% filter(time == tmax))$v )),
    as.numeric(c( (ma %>% filter(time == tmax))$v,(mp %>% filter(time == tmax))$v)),
    params$mut.strength[1],
    params$Temp,
    rep(as.character(params$web),each=(A+P)),
    c(params$degree),
    rep(nestedness_NODF(params$matrix), each=((A+P)) ),
    rep(Connectance(params$matrix), each=((A+P)) ),
    rep( (A+P),each=((A+P))),
    rep(as.character(params$var),each=(A+P))))
  
  
  colnames(species_level_data)<-c("Species","density","trait_values", "mutualism_strength", 
                                  "Temperature","webname",
                                  "Degree","Nestedness", "Connectance","Network_size",
                                  "Individual_variation")
  
  
  final_output<- list(output_network_data=output_network_data, species_level_data=species_level_data, sol=sol)
  return(final_output)
  
}


## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
## used to produce figure 1.


functional_diversity <- function(n, m, s, q = 2, delta = 0.001) {
  limits <- range(c(m - 4 * s, m + 4 * s)) # limits of integration
  taxis <- seq(limits[1], limits[2], by = delta) # trait axis
  p <- n / sum(n) # convert abundances to relative abundances
  dens <- rep(0, length(taxis)) # vector to store diversity data per grid cell
  for (i in 1:length(n)) dens <- dens + p[i]*dnorm(x=taxis, mean=m[i], sd=s[i])
  sum(dens^q)^(1/(1-q))/length(dens) # normalized functional diversity of order q
}

## code adapted from Barabas and D'Andrea 2016 Eco.Letts paper
plot_density<- function(dat) {
  dat %>%
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    theme(legend.position="none") + facet_wrap(.~type,scales = "free") %>%
    return
}


#lay out function for multiple plots
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 



#multiplot of ggplot2 figures with a common shared legend. Code taken from :https://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol =ncol , nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# plots the density distribution  at a particular timepoint. This function was used to produce figure 1.
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# moment: mean
# limits: limits of the mean trait axis which in the study are -1,1

plot_snapshot <- function(Na, Np, m, sigma, moment=0, limits=c(-1, 1), res=1001) {
  Sa <- length(Na) ## number of species
  Sp <- length(Np)
  ma<- m[1:(Sa)]
  mp<- m[Sa+1:Sp]
  sigma_a <-sigma[1:(Sa)]
  sigma_p <- sigma[Sa+1:Sp]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:Sa, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=Sa+1:Sp, trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:Sa) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:Sp) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==(Sa+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==(Sa+i))], 
                                                                mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d()+
    facet_wrap(.~species_group, nrow = 2)+
    theme(legend.title = element_text(size = 14, face = "bold"), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14, face ="bold"))+
    #geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
    #         colour="darkred", alpha=0.5, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    theme(legend.position="none") %>%
    return }




#computes the raw NODF taken from Song et al 2017 J. Animal Ecology
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}



# measures connectance of a web network
# web: interaction network
Connectance<-function(web){
  return(sum(web)/(ncol(web)*nrow(web)))}



# function for sampling competitive coefficients from random uniform distribution 
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength

mat.comp<-function(matrix){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)/Aspecies                 # scaled by number of competitors within a guild
  diag(Amatrix)<-1                                                                                           # intraspecific competition for animals
  Pmatrix<-matrix(runif(Plantspecies^2, 0.0001, 0.001), nrow=Plantspecies, ncol = Plantspecies)/Plantspecies # scaled by number of competitors within a guild
  diag(Pmatrix)<-1                                                                                           # intraspecific competion for plants
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}

# function for calculating extinction time points of every species --> stored in fact
# could be optimized by only storing time_point once for every species
calculate_extinction_tps <- function(out, network_size, fact, r) {
  all_extinct <- FALSE
  time_point <- 1
  
  while(!all_extinct) {
    species_abundances <- c((out$Na)[time_point,], (out$Np)[time_point,])
    surviving_species <- which(species_abundances > 1e-6)
    
    for(i in (surviving_species)) {
      fact$extinction_tps[[r]][i] <- time_point
    }
 
    all_extinct <- length(surviving_species) == 0

    time_point <- time_point + 1
    if(time_point > 5000) {
      break
    }
  }

  return(fact)
}




#------------------------------------------------------------------------------
#
#         Misc
#
#------------------------------------------------------------------------------

# calculate_extinction <- function(out, network_size, fact) {
#   min(which(out$Na < 1e-6))
# }

##test

##

# calculate_extinction <- function(out, network_size, fact) {
#   extinction_times_sum <- 0
#   first_extinct <- FALSE
#   all_extinct <- FALSE
#   new_extinction <- FALSE
#   species_richness <- network_size # initial richness
#   species_richness_current <- network_size
#   
#   time_point = 1
#   while(!all_extinct) {
#     species_abundances = c((out$Na)[time_point,], (out$Np)[time_point,])
#     species_richness_current <- length(which(species_abundances > 1e-6))
#     abundance_diff = species_richness_current - species_richness
#     new_extinction <- abundance_diff < 0
#     
#     if(new_extinction) {
#       extinction_times_sum = extinction_times_sum - time_point * abundance_diff # difference because abundance_diff < 0
#       
#       if(!first_extinct) {
#         fact[r,"tp_first_extinction"] = time_point
#         first_extinct <- TRUE
#       }
#       
#       species_richness <- species_richness_current
#       new_extinction <- FALSE
#       all_extinct <- species_richness == 0
#       
#     }
#     
#     time_point = time_point + 1
#     if(time_point > 5000) {
#       abundance_diff = -species_richness_current
#       extinction_times_sum = extinction_times_sum - 5000 * abundance_diff
#       break
#     }
#     #print(time_point)
#   }
#   
#   fact[r,"tp_last_extinction"] = time_point - 1
#   fact[r,"average_time_to_extinction"] = extinction_times_sum / network_size
#   
#   return(fact)
# 
