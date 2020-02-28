options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/mgh_prnp_fluid_biomarkers/')
}
library(sqldf)
library(beeswarm)

### FUNCTIONS

expand_range = function(x,by=0.5) {
  xmin = min(x,na.rm=T) - by
  xmax = max(x,na.rm=T) + by
  return (c(xmin, xmax))
}

percent = function(proportion,digits=2,format='fg') {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format=format),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

meansd = function(x, digits=1) paste(formatC(mean(x,na.rm=T),format='f',digits=digits),'±',formatC(sd(x,na.rm=T),format='f',digits=digits),sep='')

mean95 = function(x) {
  m = mean(x, na.rm=T)
  s = sd(x, na.rm=T)
  n = sum(!is.na(x))
  l95 = m - 1.96*s/sqrt(n)
  u95 = m + 1.96*s/sqrt(n)
  returnvalue = list(mean=m, l95=l95, u95=u95)
  return (returnvalue)
}

format_pval = function(p) {
  if (p > 0.0001) {
    return (formatC(p, format='fg',digits=2))
  } else {
    return (formatC(p, format='e', digits=1))
  }
}

imgmode = 'png'
# imgmode = 'pdf'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}

### CONSTANTS

participant_cats = c('-','+')
carrier_cats = c('+')
noncarrier_cats = c('-')
symptomatic_cats = c('p')
other_cohorts_color = '#C7C7C7'
lpcol = '#CD1076'

### START AN OUTPUT FILE FOR STATS TO QUOTE IN MANUSCRIPT TEXT

text_stats_path = 'figures/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T

### READ IN DATA

biomarkers = read.table('data/biomarkers.tsv', sep='\t', header=T, quote='', comment.char='', na='')
indivs = read.table('data/indivs.tsv', sep='\t', header=T, quote='', comment.char='')
params = read.table('data/params.tsv', sep='\t', header=T, quote='', comment.char='')
misc = read.table('data/misc.tsv', sep='\t', header=T, quote='', comment.char='')
prev_tr = read.table('data/prev_tr.tsv', sep='\t', header=T, quote='', comment.char='')
postlp = read.table('data/postlp_survey.tsv', sep='\t', header=T, quote='', comment.char='')

master = sqldf("
select   b.*, i.gt, p.disp, p.color, p.cat, p.x, p.panel
from     biomarkers b, indivs i, params p
where    b.deid = i.deid and i.gt = p.gt
;")

### CHECK NUMBERS

sqldf("select cat, n_visits, count(*) n from (select cat, deid, count(*) n_visits from master group by 1, 2) group by 1, 2 order by 1,2;")
sqldf("select cat, deid, count(*) n_visits from master where deid < 44 group by 1, 2 order by 1 asc, 3 desc;")

sqldf("select n_visits, count(*) n from (select deid, count(*) n_visits from biomarkers where visit > 0 group by 1) group by 1 order by 1,2;")
sqldf("select max_visit, count(*) n from (select deid, max(visit) max_visit from biomarkers where visit > 0 group by 1) group by 1 order by 1,2;")

### BEGIN FIGURE 1

imgsave(paste('figures/figure-1.',imgmode,sep=''),width=6.5*resx,height=6*resx,res=resx)

layout_matrix = as.matrix(c(1,2), nrow=2, byrow=T)
layout(layout_matrix)

vis12ppl = sqldf("
                 select   deid
                 from     master
                 where    deid in (select deid from master where visit = 1 and csf_prp is not null)
                 and      deid in (select deid from master where visit = 2 and csf_prp is not null)
                 group by 1
                 order by 1
                 ;")

vis3ppl = sqldf("
                select   deid
                from     master
                where    deid in (select deid from master where visit = 3 and csf_prp is not null)
                group by 1
                order by 1
                ;")

elisa_tr = subset(master, deid %in% vis12ppl$deid & visit %in% 1:2)
elisa_tr$l95 = elisa_tr$csf_prp - 1.96*elisa_tr$csf_prp_se
elisa_tr$u95 = elisa_tr$csf_prp + 1.96*elisa_tr$csf_prp_se

cv1 = sqldf("
            select   deid, x, gt, cat, avg(csf_prp) mean, stdev(csf_prp) sd, stdev(csf_prp)/avg(csf_prp) cv, count(*) n
            from     elisa_tr
            group by 1, 2, 3
            order by 1, 2, 3
            ;")

cv2 = sqldf("
            select   cat, avg(x) xmean, min(x) xcatmin, max(x) xcatmax, avg(cv) mean_cv, count(*) n
            from     cv1
            where    gt is not null
            group by 1
            order by 2
            ;")
cv2$mean_cv_disp = percent(cv2$mean_cv, format='f', digits=1)
# CV of individuals with anomalies:
anom_mean = mean(cv1$cv[(cv1$deid %in% misc$deid[misc$info=='csf_anomalies'])])
anom_len = length(cv1$cv[(cv1$deid %in% misc$deid[misc$info=='csf_anomalies'])])

# CV of individuals without anomalies:
no_anom_mean = mean(cv1$cv[!(cv1$deid %in% misc$deid[misc$info=='csf_anomalies'])])
no_anom_len = length(cv1$cv[!(cv1$deid %in% misc$deid[misc$info=='csf_anomalies'])])

write(paste('Figure 1 results: among ',anom_len,' individuals with anomalies noted, the mean CV was ',percent(anom_mean,digits=1,format='f'),'\n',sep=''),text_stats_path,append=T)
write(paste('Figure 1 results: among ',no_anom_len,' individuals without anomalies noted, the mean CV was ',percent(no_anom_mean,digits=1,format='f'),'\n',sep=''),text_stats_path,append=T)

betweens = seq(0.5, 10.5, by=2)
vislabs = sqldf("select (x-1)*2+visit xvis, visit from elisa_tr where x is not null and visit < 3 group by 1, 2 order by 1, 2;")
this_params = subset(params, cat != 'p')
llq = 10

par(mar=c(5,5,4,2))
plot(NA, NA, xlim=c(range(betweens)), ylim=c(0,700), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
axis(side=1, at=range(betweens), labels=NA, lwd=1, lwd.ticks=0)
mtext(side=1, line=0.25, at=vislabs$xvis, text=vislabs$visit, cex=0.9, font=3)
par(xpd=T)
mtext(side=1, line=0.25, at=0, text='visit', cex=0.9, font=3)
par(xpd=F)
mtext(side=1, line=1.5, at=this_params$x*2-.5, text=this_params$disp, font=2, col=this_params$color)
abline(v=betweens, lwd=0.5)
axis(side=2, at=(0:7)*100, lwd=1, lwd.ticks=1, las=2)
axis(side=4, at=(0:7)*100, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=2, line=3.5, text='CSF [PrP] (ng/mL)')
abline(h=llq, lwd=2, lty=3, col='#000000')
mtext(side=4, at=llq, line=0.25, las=2, text='LLQ')

axis(side=1, line=3.0, at=expand_range((elisa_tr$x[elisa_tr$cat=='-'])*2-.5, by=.75), labels=NA, tck=+0.025)
mtext(side=1, line=3.25, at=1.5, text=paste('mean CV = ',cv2$mean_cv_disp[cv2$cat=='-']))
axis(side=1, line=3.0, at=expand_range((elisa_tr$x[elisa_tr$cat=='+'])*2-.5, by=.75), labels=NA, tck=+0.025)
mtext(side=1, line=3.25, at=6.5, text=paste('mean CV = ',cv2$mean_cv_disp[cv2$cat=='+']))

for (indiv_id in unique(elisa_tr$deid)) {
  subs = subset(elisa_tr, deid==indiv_id)
  points(x=(subs$x-1)*2+subs$visit, y=subs$csf_prp, type='l', lwd=3, col=subs$color)
  points(x=(subs$x-1)*2+subs$visit, y=subs$csf_prp, pch=19, col=subs$color)
  # a few samples have two replicates identical so SEM == 0 so u95 == l95, suppress annoying zero-length arrow warnings
  suppressWarnings(arrows(x0=(subs$x-1)*2+subs$visit, x1=(subs$x-1)*2+subs$visit, y0=subs$l95, y1=subs$u95, lwd=1.5, col=subs$color, angle=90, length=0.03, code=3))
}
mtext('A', side=3, cex=2, adj = -0.1, line = 1.5)

elisa_cs = sqldf("
                 select   deid, gt, avg(csf_prp) mean_prp
                 from     elisa_tr
                 group by 1, 2
                 order by 1, 2
                 ;")

cs_aov = aov(mean_prp ~ gt, data=elisa_cs)
cs_aov_pval = as.numeric(summary(cs_aov)[[1]]['gt','Pr(>F)'])

write(paste('Figure 1 results: CSF PrP concentration showed a significant difference between genotypes (P = ',format_pval(cs_aov_pval),', one-way ANOVA) driven by lower PrP in D178N mutation carriers.\n',sep=''),text_stats_path,append=T)

summary(aov(mean_prp ~ gt, data=subset(elisa_cs, gt != 'D178N'))) # check that there is no difference without D178N
summary(aov(mean_prp ~ gt=='D178N', data=elisa_cs)) # and there is a strong difference D178N vs. all other

xmax = 365*2
yats = c((1:9)/10, 1:9, (1:9)*10)

par(mar=c(4,5,4,1))
plot(NA, NA, xlim=c(0,xmax), ylim=c(.25,4), log='y', xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=seq(0, xmax*1.1, 30.44), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0,365/2,365,365*3/2,730), labels=c(0,6,12,18,24), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=yats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.015)
axis(side=2, at=2^(-4:4), labels=percent(2^(-4:4)), lwd=0, lwd.ticks=1, tck=-0.03, las=2)
abline(h=1, lwd=1.5, lty=3)
mtext(side=1, line=2.5, text='months from first LP')
mtext(side=2, line=3.5, text='CSF [PrP] relative to 1st LP')

for (indiv in unique(prev_tr$deid)) {
  subs = subset(prev_tr, deid == indiv)
  points(subs$dayno, subs$prp_rel, col=other_cohorts_color, type='l', lwd=0.5)
}
text(x=4*30.44, y=2.5, pos=4, col=other_cohorts_color, labels='historical cohorts', font=2)

# longitudinal test-retest
tr3 = sqldf("
select   m.deid, m.visit, m.dayno_shift, m.csf_prp, m.csf_prp/v1.csf_prp prp_rel, m.gt, m.cat, m.color
from     master m, master v1
where    m.deid = v1.deid
and      v1.visit = 1
and      m.deid in (select deid from master where visit = 3 and csf_prp is not null)
order by 1, 2
;")

for (subject in unique(tr3$deid)) {
  subs = subset(tr3, deid == subject )
  points(x=subs$dayno_shift, y=subs$prp_rel, type='l', lwd=2.5, col=subs$color)
}

# compute mean cv only in the 3rd-visit people
cv1 = sqldf("
            select   deid, avg(csf_prp) mean, stdev(csf_prp) sd, count(*) n
            from     tr3
            group by 1
            order by 1
            ;")
cv2 = sqldf("
            select   avg(sd/mean) mean_cv, count(*) n
            from     cv1
            ;")
cv2

visit3_time_range = c(10,20) # - based on *unshifted* dates. otherwise: # c(floor(min(visit_days$dayno[visit_days$visit==3])/30.44), ceiling(max(visit_days$dayno[visit_days$visit==3])/30.44))

write(paste('Figure 1 results: among ',cv2$n,' individuals with longitudinal visits over ',paste(visit3_time_range,collapse='-'),' months, mean CV = ',percent(cv2$mean_cv,digits=1,format='f'),'\n',sep=''),text_stats_path,append=T)

par(xpd=T)
legend(x=18*30.44, y=8, legend=params$disp[params$cat %in% participant_cats],col=params$color[params$cat %in% participant_cats],text.col=params$color[params$cat %in% participant_cats],text.font=2,lwd=4,bty='n',cex=0.8)
par(xpd=F)

mtext('B', side=3, cex=2, adj = -0.1, line = 1.5)

dev.off() ### END FIGURE 1








### BEGIN FIGURE 2
imgsave(paste('figures/figure-2.',imgmode,sep=''),width=6.5*resx,height=6*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4),nrow=2,byrow=T)
layout(layout_matrix)

xlims = expand_range(params$x, by = 0.5)
tau_ylims = c(0,15000)
tau_short_ylims = c(0,500)
tau_yats = (0:3)*5000
tau_short_yats = (0:6)*100
tau_ylabs = tau_yats/1000 # convert pg/mL to ng/mL
tau_short_ylabs = tau_short_yats / 1000

nfl_ylims = c(0,52000)
nfl_short_ylims = c(0,2000)
nfl_yats = (0:6)*10000
nfl_short_yats = (0:4)*500
nfl_ylabs = nfl_yats / 1000 # convert pg/mL to ng/mL
nfl_short_ylabs = nfl_short_yats / 1000

par(mar=c(6,4,3,3))
plot(NA, NA, xlim=xlims, ylim=tau_ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

tau_llq = 40.6 #kitinfo3$conc[kitinfo3$cal=='CAL6']
tau_ulq = 14772.5 #kitinfo7$conc[kitinfo7$cal=='CAL1']*25/4 # poscons run at 1:25 std curve at 1:4

tau = master
# do l95 and u95 in advance so you don't get weird artifacts
tau$l95 = pmax(tau$csf_tau - 1.96*tau$csf_tau_se, tau_llq)
tau$u95 = pmin(tau$csf_tau + 1.96*tau$csf_tau_se, tau_ulq)

# tau at only most recent visit for panel A
tau_mv = sqldf("
               select   tau.*
               from     tau, (select deid, max(visit) maxvisit from tau where csf_tau is not null group by 1 order by 1) taumaxvisit
               where    tau.deid = taumaxvisit.deid
               and      tau.visit = taumaxvisit.maxvisit
               ;")

write(paste('Figure 2A (T-tau) legend: N = ',nrow(tau_mv),' samples: N = ',sum(tau_mv$cat %in% symptomatic_cats),' poscons, N = ',sum(tau_mv$cat %in% participant_cats),' participants\n',sep=''),text_stats_path,append=T)

tau_mv$xbee = 0.0

set.seed(1)
for (x in unique(tau_mv$x)) {
  rows = tau_mv$x == x & !is.na(tau_mv$csf_tau)
  tau_mv$xbee[rows] = tau_mv$x[rows] + beeswarm(csf_tau ~ x, data=tau_mv[rows,], do.plot=F, spacing=0.75, method='hex', corral='random', corralWidth=0.75)$x - 1
}

# manually fix xbee for high poscon samples for which 95%CIs are large in absolute terms
to_fix = tau_mv$pgml_av > 1000 & tau_mv$u95 - tau_mv$l95 > 200 & tau_mv$csf_tau < 6000
set.seed(2)
tau_mv$xbee[to_fix] = tau_mv$xbee[to_fix] + runif(min=-.33,max=.33,n=sum(to_fix))
to_fix_more = tau_mv$csf_tau >= 6000
tau_mv$xbee[to_fix_more] = tau_mv$xbee[to_fix_more] + runif(min=-.33,max=.33,n=sum(to_fix_more))

betweens = params$x - 0.5

limit_col = '#000000'
control_col = '#404288'
line_col = '#777777'

axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
par(xpd=T)
text(x=params$x, y=rep(0,nrow(params)), srt=45, adj=c(1,1), col=params$color, font=2, labels=params$disp)
par(xpd=F)
abline(v=betweens)
axis(side=2, at=tau_yats, labels=tau_ylabs, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='CSF T-tau (ng/mL)')
abline(h=c(tau_llq,tau_ulq), lwd=1.5, lty=2, col=limit_col)
mtext(side=4, at=tau_llq, line=0.25, las=2, col=limit_col, text='LLQ')
mtext(side=4, at=tau_ulq, line=0.25, las=2, col=limit_col, text='ULQ')
axis(side=4, at=tau_ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)

segments(x0=tau_mv$xbee, x1=tau_mv$xbee, y0=tau_mv$l95, y1=tau_mv$u95, lwd=1.5, col=tau_mv$color)
points(tau_mv$xbee, tau_mv$csf_tau, pch=20, cex=.75, col=tau_mv$color)

mtext('A', side=3, cex=2, adj = 0.0, line = 0.5)

# did CSF T-tau differ in carriers vs. non-carriers, and in participants vs. poscons?
csf_tau_ks_carrier_non = ks.test(tau_mv$csf_tau[tau_mv$cat %in% noncarrier_cats],  tau_mv$csf_tau[tau_mv$cat %in% carrier_cats], alternative='two.sided')
csf_tau_ks_symptomatic_non = ks.test(tau_mv$csf_tau[tau_mv$cat %in% participant_cats], tau_mv$csf_tau[tau_mv$cat %in% symptomatic_cats], alternative='two.sided')

participant_csftau_text = meansd(tau_mv$csf_tau[tau_mv$cat %in% participant_cats])
symptomatic_csftau_text = meansd(tau_mv$csf_tau[tau_mv$cat %in% symptomatic_cats])

write(paste('Figure 2 results: CSF T-tau, symptomatic poscons vs. all study participants: ',symptomatic_csftau_text,' vs. ',participant_csftau_text,' pg/mL, P = ',format_pval(csf_tau_ks_symptomatic_non$p.value),', two-sided Kolmogorov-Smirnov test\n',sep=''),text_stats_path,append=T)

noncarrier_csftau_text = meansd(tau_mv$csf_tau[tau_mv$cat %in% noncarrier_cats])
carrier_csftau_text = meansd(tau_mv$csf_tau[tau_mv$cat %in% carrier_cats])

write(paste('Figure 2 results: CSF T-tau, mutation non-carriers vs. carriers: ',noncarrier_csftau_text,' vs. ',carrier_csftau_text,' pg/mL, P = ',format_pval(csf_tau_ks_carrier_non$p.value),', two-sided Kolmogorov-Smirnov test\n',sep=''),text_stats_path,append=T)

# how many study visits had CSF T-tau measured?
sum(tau_mv$cat %in% participant_cats)

xlims = expand_range(params$x, by=0.5)

par(mar=c(6,4,3,3))
plot(NA, NA, xlim=xlims, ylim=nfl_ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

# nfl at only most recent visit
nfl_mv = sqldf("
               select   nfl.*
               from     master nfl, (select deid, max(visit) maxvisit from master where csf_nfl is not null group by 1 order by 1) nflmaxvisit
               where    nfl.deid = nflmaxvisit.deid
               and      nfl.visit = nflmaxvisit.maxvisit
               ;")
nfl_mv$l95 = nfl_mv$csf_nfl - 1.96*nfl_mv$csf_nfl_se
nfl_mv$u95 = nfl_mv$csf_nfl + 1.96*nfl_mv$csf_nfl_se

write(paste('Figure 2B (NfL) legend: N = ',nrow(nfl_mv),' samples: N = ',sum(nfl_mv$cat=='p'),' poscons, N = ',sum(nfl_mv$cat!='p'),' participants\n',sep=''),text_stats_path,append=T)

nfl_mv$xbee = 0.0

set.seed(1)
for (x in unique(nfl_mv$x)) {
  rows = nfl_mv$x == x
  nfl_mv$xbee[rows] = nfl_mv$x[rows] + beeswarm(csf_nfl ~ x, data=nfl_mv[rows,], do.plot=F, spacing=0.75, method='hex', corral='random', corralWidth=0.75)$x - 1
}

betweens = params$x - 0.5

nfl_llq = 100 * 2
nfl_ulq = 10000 * 5

axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
par(xpd=T)
text(x=params$x, y=rep(0,nrow(params)), srt=45, adj=c(1,1), col=params$color, font=2, labels=params$disp)
par(xpd=F)

abline(v=betweens)
axis(side=2, at=nfl_yats, labels=nfl_ylabs, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='CSF NfL (ng/mL)')
abline(h=c(nfl_llq,nfl_ulq), lwd=1.5, lty=2, col=limit_col)
mtext(side=4, at=nfl_llq, line=0.25, las=2, col=limit_col, text='LLQ')
mtext(side=4, at=nfl_ulq, line=0.25, las=2, col=limit_col, text='ULQ')
axis(side=4, at=nfl_ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)

segments(x0=nfl_mv$xbee, x1=nfl_mv$xbee, y0=nfl_mv$l95, y1=nfl_mv$u95, lwd=1.5, col=nfl_mv$color)
points(nfl_mv$xbee, nfl_mv$csf_nfl, pch=20, cex=.75, col=nfl_mv$color)

mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

# did CSF NfL differ in carriers vs. non-carriers, and in participants vs. poscons?
csf_nfl_ks_carrier_non = ks.test(nfl_mv$csf_nfl[nfl_mv$cat %in% noncarrier_cats],  nfl_mv$csf_nfl[nfl_mv$cat %in% carrier_cats], alternative='two.sided')
csf_nfl_ks_symptomatic_non = ks.test(nfl_mv$csf_nfl[nfl_mv$cat %in% participant_cats], nfl_mv$csf_nfl[nfl_mv$cat %in% symptomatic_cats], alternative='two.sided')

participant_csfnfl_text = meansd(nfl_mv$csf_nfl[nfl_mv$cat %in% participant_cats])
symptomatic_csfnfl_text = meansd(nfl_mv$csf_nfl[nfl_mv$cat %in% symptomatic_cats])

write(paste('Figure 2 results: CSF NfL, symptomatic poscons vs. all study participants: ',symptomatic_csfnfl_text,' vs. ',participant_csfnfl_text,' pg/mL, P = ',format_pval(csf_nfl_ks_symptomatic_non$p.value),', two-sided Kolmogorov-Smirnov test\n',sep=''),text_stats_path,append=T)

noncarrier_csfnfl_text = meansd(nfl_mv$csf_nfl[nfl_mv$cat %in% noncarrier_cats])
carrier_csfnfl_text = meansd(nfl_mv$csf_nfl[nfl_mv$cat %in% carrier_cats])

write(paste('Figure 2 results: CSF NfL, mutation non-carriers vs. carriers: ',noncarrier_csfnfl_text,' vs. ',carrier_csfnfl_text,' pg/mL, P = ',format_pval(csf_nfl_ks_carrier_non$p.value),', two-sided Kolmogorov-Smirnov test\n',sep=''),text_stats_path,append=T)

sum(nfl_mv$cat %in% participant_cats)

### panels C & D - longitudinal tau & NfL where avails

tau_longit = sqldf("
                   select   *
                   from     master t
                   where    t.deid in (select deid from master where visit = 3 and csf_tau is not null)
                   and      csf_tau is not null
                   order by deid, visit
                   ;")
m = lm(csf_tau ~ deid + dayno_shift, data=tau_longit)
tau_longit_pval = summary(m)$coefficients['dayno_shift','Pr(>|t|)']
write(paste('Figure 2C results: CSF T-tau longitudinal P value: ',format_pval(tau_longit_pval),', linear regression, N = ',length(unique(tau_longit$deid)),' individuals\n',sep=''),text_stats_path,append=T)

tau_longit_cv_calc = sqldf("
                           select   deid, avg(csf_tau) mean_tau, stdev(csf_tau) sd_tau, count(*) n, stdev(csf_tau)/avg(csf_tau) cv_tau
                           from     tau_longit
                           group by 1
                           order by 1
                           ;")
tau_longit_cv = mean(tau_longit_cv_calc$cv)
write(paste('Figure 2C results: CSF T-tau longitudinal CV: ',percent(tau_longit_cv, digits=1, format='f'),', N = ',nrow(tau_longit_cv_calc),' individuals.\n',sep=''),text_stats_path,append=T)

m = lm(csf_tau ~ deid + dayno_shift, data=subset(tau_longit, deid != '17'))
tau_longit_pval = summary(m)$coefficients['dayno_shift','Pr(>|t|)']
write(paste('Figure 2C results: CSF T-tau 3-visit longitudinal P value without the one mutation-negative individual who showed a suggestive upward trend: ',format_pval(tau_longit_pval),', linear regression\n',sep=''),text_stats_path,append=T)

par(mar=c(4,4,3,3))
xmax = 365*2
plot(NA, NA, xlim=c(0,xmax), ylim=tau_short_ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=seq(0, xmax*1.1, 30.44), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0,365/2,365,365*3/2,730), labels=c(0,6,12,18,24), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=tau_short_yats, labels=tau_short_ylabs, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='CSF T-tau (ng/mL)')
abline(h=c(tau_llq,tau_ulq), lwd=1.5, lty=2, col=limit_col)
mtext(side=4, at=tau_llq, line=0.25, las=2, col=limit_col, text='LLQ')
mtext(side=4, at=tau_ulq, line=0.25, las=2, col=limit_col, text='ULQ')
axis(side=4, at=tau_short_ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=1, line=2.5, text='months from first LP')

for (each_indiv in unique(tau_longit$deid)) {
  subs = subset(tau_longit, deid==each_indiv)
  points(subs$dayno_shift, subs$csf_tau, type='l', lwd=3, col=subs$color)
}

mtext('C', side=3, cex=2, adj = -0.1, line = 1.5)

nfl_longit = sqldf("
                   select   *
                   from     master n
                   where    n.deid in (select deid from master where visit = 3 and csf_nfl is not null)
                   and      csf_nfl is not null
                   order by deid, visit
                   ;")

m = lm(csf_nfl ~ deid + dayno_shift, data=nfl_longit)
nfl_longit_pval = summary(m)$coefficients['dayno_shift','Pr(>|t|)']
write(paste('Figure 2D results: CSF NfL longitudinal P value: ',format_pval(nfl_longit_pval),', linear regression, N = ',length(unique(nfl_longit$deid)),' individuals.\n',sep=''),text_stats_path,append=T)

nfl_longit_cv_calc = sqldf("
                           select   deid, avg(csf_nfl) mean_nfl, stdev(csf_nfl) sd_nfl, count(*) n, stdev(csf_nfl)/avg(csf_nfl) cv
                           from     nfl_longit
                           group by 1
                           order by 1
                           ;")
nfl_longit_cv = mean(nfl_longit_cv_calc$cv)
write(paste('Figure 2D results: CSF NfL longitudinal CV: ',percent(nfl_longit_cv, digits=1, format='f'),', N = ',nrow(nfl_longit_cv_calc),' individuals.\n',sep=''),text_stats_path,append=T)

par(mar=c(4,4,3,3))
xmax = 365*2
plot(NA, NA, xlim=c(0,xmax), ylim=nfl_short_ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=seq(0, xmax*1.1, 30.44), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0,365/2,365,365*3/2,730), labels=c(0,6,12,18,24), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=nfl_short_yats, labels=nfl_short_ylabs, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='CSF NfL (ng/mL)')
abline(h=c(nfl_llq,nfl_ulq), lwd=1.5, lty=2, col=limit_col)
mtext(side=4, at=nfl_llq, line=0.25, las=2, col=limit_col, text='LLQ')
mtext(side=4, at=nfl_ulq, line=0.25, las=2, col=limit_col, text='ULQ')
axis(side=4, at=nfl_short_ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=1, line=2.5, text='months from first LP')

for (each_indiv in unique(nfl_longit$deid)) {
  subs = subset(nfl_longit, deid==each_indiv)
  points(subs$dayno_shift, subs$csf_nfl, type='l', lwd=3, col=subs$color)
}

mtext('D', side=3, cex=2, adj = -0.1, line = 1.5)

dev.off() ### END FIGURE 2








### BEGIN FIGURE 3
imgsave(paste('figures/figure-3.',imgmode,sep=''),width=6.5*resx,height=3.25*resx,res=resx)

layout_matrix = matrix(c(1,2),nrow=1,byrow=T)
layout(layout_matrix)

nfl_ref_ranges = data.frame(value=c(222.8,204.8,296.6,150),
                            study=c('Steinacker 2016', 'Kovacs 2017', 'Thompson 2018', 'Staffaroni 2019'))

plasma = master[!is.na(master$plasma_nfl) | !is.na(master$plasma_tau),c('deid','visit','cat','disp','x','color','dayno_shift','plasma_nfl','plasma_tau')]

# most recent visit only
plasma_indiv = sqldf("
                     select   p.deid, p.visit, p.cat, p.disp, p.x, p.color, p.plasma_tau tau, p.plasma_nfl nfl
                     from     plasma p, (select deid, max(visit) maxvisit from plasma where plasma_nfl is not null group by 1 order by 1) l
                     where    p.deid = l.deid and p.visit = l.maxvisit
                     group by 1, 2, 3, 4, 5
                     order by 1, 2
                     ;")
plasma_nfl_llq = 2.7 # pg/mL
plasma_indiv$nfl[plasma_indiv$nfl < plasma_nfl_llq] = plasma_nfl_llq

xlims = expand_range(params$x[params$cat %in% participant_cats],.5)
plasma_tau_ylims = c(0,30)

plasma_nfl_ylims = c(0,300)
par(mar=c(4,4,3,3))
plot(NA, NA, xlim=xlims, ylim=plasma_nfl_ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

plasma_indiv$xbee = plasma_indiv$x # set a default (so the NA rows do not get lost)
set.seed(1)
for (x in unique(plasma_indiv$x)) {
  rows = plasma_indiv$x == x & !is.na(plasma_indiv$tau)
  plasma_indiv$xbee[rows] = plasma_indiv$x[rows] + beeswarm(tau ~ x, data=plasma_indiv[rows,], do.plot=F, spacing=0.75, method='hex', corral='random', corralWidth=0.75)$x - 1
}

betweens = params$x - 0.5
xlabs = sqldf("select disp, color, avg(x) xmid from plasma_indiv group by 1, 2 order by 1;")

# did plasma tau differ in carriers vs. non-carriers
plasma_tau_ks_carrier_non     = ks.test(plasma_indiv$tau[plasma_indiv$cat %in% noncarrier_cats],  plasma_indiv$tau[plasma_indiv$cat %in% carrier_cats], alternative='two.sided')

noncarrier_plasmatau_text = meansd(plasma_indiv$tau[plasma_indiv$cat %in% noncarrier_cats])
carrier_plasmatau_text    = meansd(plasma_indiv$tau[plasma_indiv$cat %in% carrier_cats])

write(paste('Figure 3 results: plasma tau, mutation non-carriers vs. carriers: ',noncarrier_plasmatau_text,' vs. ',carrier_plasmatau_text,' pg/mL, P = ',format_pval(plasma_tau_ks_carrier_non$p.value),', two-sided Kolmogorov-Smirnov test\n',sep=''),text_stats_path,append=T)

axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
par(xpd=T)
text(x=xlabs$xmid, y=rep(0,nrow(xlabs)), srt=45, adj=c(1,1), col=xlabs$color, font=2, labels=xlabs$disp)
par(xpd=F)
abline(v=betweens)
abline(h=nfl_ref_ranges$value, lwd=0.5, lty=3)
par(xpd=T)
text(x=6,y=mean(nfl_ref_ranges$value),label='symptomatic\nreference ranges',srt=270)
par(xpd=F)
axis(side=2, at=(0:6)*50, labels=(0:6)*50, lwd=1, lwd.ticks=1, las=2)
abline(h=plasma_nfl_llq, lwd=0.5, lty=2)
mtext(side=4,line=0.5,at=plasma_nfl_llq,las=2,text='LLQ')
mtext(side=2, line=2.5, text='plasma NfL (pg/mL)')
points(plasma_indiv$xbee, plasma_indiv$nfl, pch=20, cex=.75, col=plasma_indiv$color)
mtext('A', side=3, cex=2, adj = 0.0, line = 0.5)

write(paste('Figure 3A legend: N= ',sum(!is.na(plasma_indiv$nfl)),' individuals.\n',sep=''),text_stats_path,append=T)

# did plasma NfL differ in carriers vs. non-carriers
plasma_nfl_ks_carrier_non     = ks.test(plasma_indiv$nfl[plasma_indiv$cat %in% noncarrier_cats],  plasma_indiv$nfl[plasma_indiv$cat %in% carrier_cats], alternative='two.sided')

noncarrier_plasmanfl_text = meansd(plasma_indiv$nfl[plasma_indiv$cat %in% noncarrier_cats])
carrier_plasmanfl_text    = meansd(plasma_indiv$nfl[plasma_indiv$cat %in% carrier_cats])

write(paste('Figure 3 results: plasma NfL, mutation non-carriers vs. carriers: ',noncarrier_plasmanfl_text,' vs. ',carrier_plasmanfl_text,' pg/mL, P = ',format_pval(plasma_nfl_ks_carrier_non$p.value),', two-sided Kolmogorov-Smirnov test\n',sep=''),text_stats_path,append=T)

plasma_longit_3visit = sqldf("select * from plasma where deid in (select deid from plasma where visit = 3);")

par(mar=c(4,4,3,3))
xmax = 365*2
ymax = 20
plot(NA, NA, xlim=c(0,xmax), ylim=c(0,ymax), xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=seq(0, xmax*1.1, 30.44), labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=c(0,365/2,365,365*3/2,730), labels=c(0,6,12,18,24), lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=(0:4)*5, labels=(0:4)*5, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.0, text='plasma NfL (pg/mL)')
abline(h=plasma_nfl_llq, lwd=0.5, lty=2)
mtext(side=4,line=0.5,at=plasma_nfl_llq,las=2,text='LLQ')
mtext(side=1, line=2.5, text='months from first visit')

for (each_indiv in unique(plasma_longit_3visit$deid)) {
  subs = subset(plasma_longit_3visit, deid==each_indiv)
  points(subs$dayno_shift, subs$plasma_nfl, type='l', lwd=3, col=subs$color)
}

m = lm(nfl ~ dayno_shift + deid, data=plasma_longit_3visit)
summary(m)
plasma_nfl_longit_pval = summary(m)$coefficients['dayno_shift','Pr(>|t|)']
write(paste('Figure 3 results: plasma NfL 3-visit longitudinal P value: ',format_pval(nfl_longit_pval),', linear regression\n',sep=''),text_stats_path,append=T)
mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

write(paste('Figure 3B legend: N= ',length(unique(plasma_longit_3visit$deid)),' individuals.\n',sep=''),text_stats_path,append=T)

dev.off() ### END FIGURE 3






### BEGIN FIGURE 4
imgsave(paste('figures/figure-4.',imgmode,sep=''),width=6.5*resx,height=7.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,5,6),nrow=3,byrow=T)
layout(layout_matrix)

rtq_xlims = c(0,24)
rtq_xats = 0:24
rtq_xbigs = seq(0,24,6)
rtq_ylims = c(-0.025,1.05)
rtq_yats = (0:4)/4

rtq_kinetics = read.table('data/rtq_sha_kinetic.tsv',sep='\t',header=T)
rtq_kinetics$cat = master$cat[match(rtq_kinetics$deid, master$deid_sample)]
rtq_kinetics$gt = master$gt[match(rtq_kinetics$deid, master$deid_sample)]

for (panel in 1:6) {
  i = which(params$panel==panel)
  
  par(mar=c(4,5,4,2))
  plot(NA, NA, xlim=rtq_xlims, ylim=rtq_ylims, xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
  axis(side=1, at=rtq_xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=rtq_xbigs, labels=rtq_xbigs, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=2, at=rtq_yats, labels=percent(rtq_yats), las=2)
  mtext(side=1, line=2.5, text='reaction time (h)')
  mtext(side=2, line=3, text='normalized fluorescence')
  
  for (each_sample in unique(rtq_kinetics$deid[rtq_kinetics$gt==params$gt[i]])) {
    subs = subset(rtq_kinetics, deid == each_sample)
    points(x=subs$hours, y=subs$fluor, type='l', lwd=3, col=params$color[i])
  }
  
  # compute summary stats for this category
  positives = sum(master$gt == params$gt[i] & master$rtq_sha_overall_call == '+' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)
  negatives = sum(master$gt == params$gt[i] & master$rtq_sha_overall_call == '-' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)
  positive_rate_message = paste(positives, positives + negatives, sep='/')
  
  mtext(side=3, line=2, text=params$disp[i], col=params$color[i], font=2)
  mtext(side=3, line=0.5, text=positive_rate_message, col='#000000', font=1)
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
}

dev.off() ### END FIGURE 4

# double check numbers
length(unique(rtq_kinetics$deid))
sum(!is.na(master$rtq_sha_overall_call))

poscon_positives = sum(master$cat == 'p' & master$rtq_sha_overall_call == '+' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)
poscon_negatives = sum(master$cat == 'p' & master$rtq_sha_overall_call == '-' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)

write(paste('Figure 4 results: SHaPrP RT-QuIC sensitivity to positive controls was: ',percent(poscon_positives/(poscon_positives+poscon_negatives)),'\n',sep=''),text_stats_path,append=T)





### BEGIN FIGURE S1
imgsave(paste('figures/figure-s1.',imgmode,sep=''),width=3.25*resx,height=3.25*resx,res=resx)

par(mar=c(3,4,2,1))
plot(NA, NA, xlim=c(0.75, 2.25), ylim=c(0,100), xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=c(0,1,2,3), labels=c('','pre-LP','post-LP',''), lwd=1, lwd.ticks=0)
axis(side=2, at=(0:4)*25, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='subjective rating')
points(x=rep(1,nrow(postlp)), y=postlp$beforelp_v1, pch=20, col=lpcol)
points(x=rep(2,nrow(postlp)), y=postlp$futurelp_v1, pch=20, col=lpcol)
segments(x0=rep(1,nrow(postlp)), x1=rep(2,nrow(postlp)), y0=postlp$beforelp_v1, y1=postlp$futurelp_v1, lwd=1.5, col=lpcol)

points(x=0.9,  y= mean95(postlp$beforelp_v1)$mean, pch=15, col='black')
arrows(x0=0.9, y0=mean95(postlp$beforelp_v1)$l95, y1=mean95(postlp$beforelp_v1)$u95, angle=90, length=0.05, code=3)
points(x=2.1,  y= mean95(postlp$futurelp_v1)$mean, pch=15, col='black')
arrows(x0=2.1, y0=mean95(postlp$futurelp_v1)$l95, y1=mean95(postlp$futurelp_v1)$u95, angle=90, length=0.05, code=3)

tobj = t.test(postlp$beforelp_v1, postlp$futurelp_v1, paired=T, alternative='two.sided')
write(paste('Figure S1 legend: N = ',nrow(postlp),', nominal P =  ',formatC(tobj$p.value, digits=2),', two-sided paired t test\n',sep=''),text_stats_path,append=T)

dev.off()












### BEGIN FIGURE S2
imgsave(paste('figures/figure-s2.',imgmode,sep=''),width=6.5*resx,height=7.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,byrow=T)
layout(layout_matrix)

panel = 1

gothcsf_indiv = sqldf("
                      select   s.deid, s.gt, s.disp, s.cat, s.x, s.color, 
                               s.gothenburg_csf_tau tau, s.gothenburg_csf_nfl nfl, s.gothenburg_csf_gfap gfap
                      from     master s, (select deid, max(visit) maxvisit from master where gothenburg_csf_tau is not null group by 1 order by 1) l
                      where    s.deid = l.deid and s.visit = l.maxvisit
                      group by 1, 2, 3, 4, 5
                      order by 1, 2
                      ;")
xlims = expand_range(params$x[params$cat %in% participant_cats],.5)

yats = (0:20)*100
ybigs = (0:4)*500
ybigs_labs = formatC(ybigs, big.mark=',')

par(mar=c(6,5,3,2))
plot(NA, NA, xlim=xlims, ylim=c(0,1000), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

gothcsf_indiv$xbee = gothcsf_indiv$xcat # set a default so NA do not get lost
set.seed(1)
for (x in unique(gothcsf_indiv$x)) {
  rows = gothcsf_indiv$x == x & !is.na(gothcsf_indiv$tau)
  gothcsf_indiv$xbee[rows] = gothcsf_indiv$x[rows] + beeswarm(tau ~ x, data=gothcsf_indiv[rows,], do.plot=F, spacing=0.75, method='hex', corral='random', corralWidth=0.75)$x - 1
}

betweens = params$x - 0.5

par(new=TRUE) # necessary b/c we already called plot() above to frame the xbee calculation

for (marker in c('tau','NfL','GFAP')) {
  
  mkr = tolower(marker)
  ymax = ceiling(max(gothcsf_indiv[,mkr],na.rm=T)*1.1/100)*100
  
  plot(NA, NA, xlim=xlims, ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd=1, lwd.ticks=0)
  xlabs = sqldf("select disp, color, avg(x) xmid from gothcsf_indiv group by 1, 2 order by 1;")
  par(xpd=T)
  text(x=xlabs$xmid, y=rep(0,nrow(xlabs)), srt=45, adj=c(1,1), col=xlabs$color, font=2, labels=xlabs$disp)
  par(xpd=F)
  abline(v=betweens)
  axis(side=2, at=yats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025, las=2)
  axis(side=2, at=ybigs, labels=ybigs_labs, lwd=1, lwd.ticks=1, tck=-0.050, las=2)
  marker_label = marker
  if (marker == 'tau') {
    marker_label = 'T-tau'
  }
  mtext(side=2, line=3.5, text=paste('CSF ',marker_label,' (pg/mL)',sep=''))
  points(gothcsf_indiv$xbee, gothcsf_indiv[,mkr], pch=20, cex=.75, col=gothcsf_indiv$color)
  mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 0.5)
  
  write(paste('Figure S2',LETTERS[panel],' legend: ',mkr,' cross-sectional comparison based on N = ',length(unique(gothcsf_indiv$deid)),' individuals.\n',sep=''),text_stats_path,append=T)
  
  panel = panel + 1
}

# 2-visit test-retest CV plots
gothcsf_tau_indiv = sqldf("select deid, avg(gothenburg_csf_tau) mean_tau, stdev(gothenburg_csf_tau) sd_tau, count(*) n from master where gothenburg_csf_tau is not null group by 1 having count(*) > 1 order by 1;")
gothcsf_tau_indiv$cv = gothcsf_tau_indiv$sd_tau / gothcsf_tau_indiv$mean_tau
gothcsf_tau_cv = mean(gothcsf_tau_indiv$cv)

gothcsf_nfl_indiv = sqldf("select deid, avg(gothenburg_csf_nfl) mean_nfl, stdev(gothenburg_csf_nfl) sd_nfl, count(*) n from master where gothenburg_csf_nfl is not null  group by 1 having count(*) > 1 order by 1;")
gothcsf_nfl_indiv$cv = gothcsf_nfl_indiv$sd_nfl / gothcsf_nfl_indiv$mean_nfl
gothcsf_nfl_cv = mean(gothcsf_nfl_indiv$cv)

gothcsf_gfap_indiv = sqldf("select deid, avg(gothenburg_csf_gfap) mean_gfap, stdev(gothenburg_csf_gfap) sd_gfap, count(*) n from master where gothenburg_csf_gfap is not null  group by 1 having count(*) > 1 order by 1;")
gothcsf_gfap_indiv$cv = gothcsf_gfap_indiv$sd_gfap / gothcsf_gfap_indiv$mean_gfap
gothcsf_gfap_cv = mean(gothcsf_gfap_indiv$cv)

par(mar=c(4,5,4,2))
for (marker in c('tau','NfL','GFAP')) {
  
  mkr = tolower(marker)
  colname = paste0('gothenburg_csf_',mkr)
  ymax = ceiling(max(master[,colname],na.rm=T)*1.1/100)*100
  
  par(mar=c(6,5,3,2))
  plot(NA, NA, xlim=c(0.5, 2.5), ylim=c(0,ymax), axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=c(0.5,2.5), labels=NA, lwd=1, lwd.ticks=0)
  axis(side=1, at=c(1,2), labels=c(1,2), lwd=0, lwd.ticks=0)
  mtext(side=1, line=3, text='study visit')
  axis(side=2, at=yats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025, las=2)
  axis(side=2, at=ybigs, labels=ybigs_labs, lwd=1, lwd.ticks=1, tck=-0.050, las=2)
  marker_label = marker
  if (marker == 'tau') {
    marker_label = 'T-tau'
  }
  mtext(side=2, line=3.5, text=paste('CSF ',marker_label,' (pg/mL)',sep=''))
  two_visit_indivs = sqldf("select deid from master where gothenburg_csf_tau is not null group by 1 having count(*) > 1 order by 1;")
  for (individ in two_visit_indivs$deid) {
    s = master[master$deid==individ,]
    points(x=s$visit, y=s[,colname], col=s$color, pch=19)
    segments(x0=1,x1=2,y0=s[s$visit==1,colname],y1=s[s$visit==2,colname], col=s$color, lwd=3)
  }
  mtext(side=3, line=0, text=paste('mean CV: ',percent(get(paste('gothcsf_',mkr,'_cv',sep=''))),sep=''),cex=0.9)
  mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 1.5)
  
  write(paste('Figure S2',LETTERS[panel],' legend: ',mkr,' test-retest based on N = ',sum(!is.na(master[master$deid %in% two_visit_indivs$deid,colname])),' samples from N = ',length(unique(two_visit_indivs$deid)),' individuals.\n',sep=''),text_stats_path,append=T)
  
  panel = panel + 1
}

nflcomp = sqldf("
                select   deid, visit, csf_nfl broad, gothenburg_csf_nfl goth, color
                from     master m
                where    csf_nfl is not null and gothenburg_csf_nfl is not null
                order by 1, 2
                ;")
range(nflcomp$broad)
range(nflcomp$goth)
nrow(nflcomp)
length(unique(nflcomp$deid))
write(paste('Figure S2G legend: CSF NfL comparison based on N = ',nrow(nflcomp),' samples from N = ',length(unique(nflcomp$deid)),' individuals.\n',sep=''),text_stats_path,append=T)

taucomp = sqldf("
                select   deid, visit, csf_tau broad, gothenburg_csf_tau goth, color
                from     master m
                where    csf_tau is not null and gothenburg_csf_tau is not null
                order by 1, 2
                ;")
range(taucomp$broad)
range(taucomp$goth)
nrow(taucomp)
length(unique(taucomp$deid))
write(paste('Figure S2H legend: CSF T-tau comparison based on N = ',nrow(taucomp),' samples from N = ',length(unique(taucomp$deid)),' individuals.\n',sep=''),text_stats_path,append=T)

nfl_lims = c(0,2200)
tau_lims = c(0,650)

ats = (0:21)*100
bigs = (0:5)*500
par(mar=c(4,5,3,2))
plot(NA, NA, xlim=nfl_lims, ylim=nfl_lims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=bigs, labels=bigs, lwd=1, lwd.ticks=1, tck=-0.050)
mtext(side=1, line=2.5, text='Broad CSF NfL pg/mL')
axis(side=2, at=ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=bigs, labels=bigs, lwd=1, lwd.ticks=1, tck=-0.050, las=2)
mtext(side=2, line=3.5, text='Gothenburg CSF NfL pg/mL')
abline(a=0, b=1)
points(nflcomp$broad, nflcomp$goth, pch=20, col=nflcomp$color)
spearman = cor.test(nflcomp$broad, nflcomp$goth, method='spearman')
msg = paste0('r = ',formatC(spearman$estimate,format='f',digits=2),', P = ',format_pval(spearman$p.value))
mtext(side=3, line=0, text=msg)
mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 1.5)
panel = panel + 1

ats = (0:21)*100
bigs = (0:5)*500
par(mar=c(4,5,3,2))
plot(NA, NA, xlim=tau_lims, ylim=tau_lims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=bigs, labels=bigs, lwd=1, lwd.ticks=1, tck=-0.050)
mtext(side=1, line=2.5, text='Broad CSF T-tau pg/mL')
axis(side=2, at=ats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=2, at=bigs, labels=bigs, lwd=1, lwd.ticks=1, tck=-0.050, las=2)
mtext(side=2, line=3.5, text='Gothenburg CSF T-tau pg/mL')
abline(a=0, b=1)
points(taucomp$broad, taucomp$goth, pch=20, col=taucomp$color)
spearman = cor.test(taucomp$broad, taucomp$goth, method='spearman')
msg = paste0('r = ',formatC(spearman$estimate,format='f',digits=2),', P = ',format_pval(spearman$p.value))
mtext(side=3, line=0, text=msg)
mtext(LETTERS[panel], side=3, cex=2, adj = 0.0, line = 1.5)
panel = panel + 1

# blank final panel
plot(NA, NA, xlim=0:1, ylim=0:1, axes=F, ann=F, xaxs='i', yaxs='i')

dev.off() ### END FIGURE S2










### BEGIN FIGURE S3
imgsave(paste('figures/figure-s3.',imgmode,sep=''),width=6.5*resx,height=3.0*resx,res=resx)

layout_matrix = matrix(c(1,1,2,2,3), nrow=1, byrow=T)
#layout_matrix = matrix(c(1,2),nrow=1,byrow=T)
layout(layout_matrix)

plasma_tau_indiv = sqldf("select deid, avg(plasma_tau) mean_tau, stdev(plasma_tau) sd_tau, count(*) n from master where plasma_tau is not null and visit in (1,2)  group by 1 having count(*) > 1 order by 1;")
plasma_tau_indiv$cv = plasma_tau_indiv$sd_tau / plasma_tau_indiv$mean_tau
plasma_tau_cv = mean(plasma_tau_indiv$cv)
plasma_tau_llq = 1.22
write(paste('Figure S3 legend: plasma T-tau short-term mean CV = ',percent(plasma_tau_cv),', N = ',nrow(plasma_tau_indiv),'\n',sep=''),text_stats_path,append=T)

plasma_longit = sqldf('select * from master where plasma_tau is not null and visit in (1,2);')
two_tau_indivs = sqldf("select deid from plasma_longit where plasma_tau is not null and visit in (1,2) group by 1 having count(*) > 1 order by 1;")

par(mar=c(5,5,4,2))
plot(NA, NA, xlim=c(0.5, 2.5), ylim=c(0,30), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=c(0.5,2.5), labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=c(1,2), labels=c(1,2), lwd=0, lwd.ticks=0)
mtext(side=1, line=3, text='study visit')
axis(side=2, at=c(0:6)*5, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='plasma tau (pg/mL)')
for (deid in two_tau_indivs$deid) {
  rows = plasma_longit$deid==deid
  points(x=plasma_longit$visit[rows], plasma_longit$plasma_tau[rows], col=plasma_longit$color[rows], pch=19)
  points(x=plasma_longit$visit[rows], plasma_longit$plasma_tau[rows], col=plasma_longit$color[rows], type='l', lwd=3)
}
mtext(side=3, line=0, text=paste('mean CV: ',percent(plasma_tau_cv),sep=''),cex=0.9)
abline(h=plasma_tau_llq, lwd=1, lty=3)
axis(side=4, at=plasma_tau_llq, labels='LLQ', lwd=0, lwd.ticks=0, las=2)
mtext('A', side=3, cex=2, adj = 0.0, line = 1.5)

plasma_nfl_indiv = sqldf("select deid, avg(plasma_nfl) mean_nfl, stdev(plasma_nfl) sd_nfl, count(*) n from master where plasma_nfl is not null and visit in (1,2) group by 1 having count(*) > 1 order by 1;")
plasma_nfl_indiv$cv = plasma_nfl_indiv$sd_nfl / plasma_nfl_indiv$mean_nfl
plasma_nfl_cv = mean(plasma_nfl_indiv$cv)
plasma_nfl_llq = 2.7
write(paste('Figure S3 legend: plasma NfL short-term mean CV = ',percent(plasma_nfl_cv),', N = ',nrow(plasma_nfl_indiv),'\n',sep=''),text_stats_path,append=T)

plasma_longit = sqldf('select * from master where plasma_nfl is not null and visit in (1,2);')
two_nfl_indivs = sqldf("select deid from plasma_longit where plasma_nfl is not null and visit in (1,2) group by 1 having count(*) > 1 order by 1;")

plasma_longit$plasma_nfl[plasma_longit$plasma_nfl < plasma_nfl_llq] = plasma_nfl_llq

par(mar=c(5,5,4,2))
plot(NA, NA, xlim=c(0.5, 2.5), ylim=c(0,30), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=c(0.5,2.5), labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=c(1,2), labels=c(1,2), lwd=0, lwd.ticks=0)
mtext(side=1, line=3, text='study visit')
axis(side=2, at=c(0:6)*5, lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='plasma NfL (pg/mL)')
for (deid in two_nfl_indivs$deid) {
  rows = plasma_longit$deid==deid
  points(x=plasma_longit$visit[rows], plasma_longit$plasma_nfl[rows], col=plasma_longit$color[rows], pch=19)
  points(x=plasma_longit$visit[rows], plasma_longit$plasma_nfl[rows], col=plasma_longit$color[rows], type='l', lwd=3)
}
mtext(side=3, line=0, text=paste('mean CV: ',percent(plasma_nfl_cv),sep=''),cex=0.9)
abline(h=nfl_ref_ranges$value, lwd=0.5, lty=3)
abline(h=plasma_nfl_llq, lwd=1, lty=3)
axis(side=4, at=plasma_nfl_llq, labels='LLQ', lwd=0, lwd.ticks=0, las=2)
mtext('B', side=3, cex=2, adj = 0.0, line = 1.5)

par(mar=c(5,1,4,1))
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', axes=F, ann=F)
legend(x=0, y=1, legend=params$disp[params$cat %in% participant_cats],col=params$color[params$cat %in% participant_cats],text.col=params$color[params$cat %in% participant_cats],text.font=2,lwd=4,bty='n',cex=1)

dev.off() ### END FIGURE S3









### BEGIN FIGURE S4
imgsave(paste('figures/figure-s4.',imgmode,sep=''),width=6.5*resx,height=7.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,5,6),nrow=3,byrow=T)
layout(layout_matrix)

rtq_xlims = c(0,24)
rtq_xats = 0:24
rtq_xbigs = seq(0,24,6)
rtq_ylims = c(-0.025,1.05)
rtq_yats = (0:4)/4

rtq_kinetics = read.table('data/rtq_bv_kinetic.tsv',sep='\t',header=T)
rtq_kinetics$cat = master$cat[match(rtq_kinetics$deid, master$deid_sample)]
rtq_kinetics$gt = master$gt[match(rtq_kinetics$deid, master$deid_sample)]

for (panel in 1:6) {
  i = which(params$panel==panel)
  
  par(mar=c(4,5,4,2))
  plot(NA, NA, xlim=rtq_xlims, ylim=rtq_ylims, xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
  axis(side=1, at=rtq_xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=rtq_xbigs, labels=rtq_xbigs, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=2, at=rtq_yats, labels=percent(rtq_yats), las=2)
  mtext(side=1, line=2.5, text='reaction time (h)')
  mtext(side=2, line=3, text='normalized fluorescence')
  
  for (each_sample in unique(rtq_kinetics$deid[rtq_kinetics$gt==params$gt[i]])) {
    subs = subset(rtq_kinetics, deid == each_sample)
    points(x=subs$hours, y=subs$fluor, type='l', lwd=3, col=params$color[i])
  }
  
  # compute summary stats for this category
  positives = sum(master$gt == params$gt[i] & master$rtq_bv_overall_call == '+' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)
  negatives = sum(master$gt == params$gt[i] & master$rtq_bv_overall_call == '-' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)
  positive_rate_message = paste(positives, positives + negatives, sep='/')
  
  mtext(side=3, line=2, text=params$disp[i], col=params$color[i], font=2)
  mtext(side=3, line=0.5, text=positive_rate_message, col='#000000', font=1)
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 1.5)
}

dev.off() ### END FIGURE S4

# double check numbers
length(unique(rtq_kinetics$deid))
sum(!is.na(master$rtq_bv_overall_call))

poscon_positives = sum(master$cat == 'p' & master$rtq_bv_overall_call == '+' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)
poscon_negatives = sum(master$cat == 'p' & master$rtq_bv_overall_call == '-' & master$deid_sample %in% rtq_kinetics$deid, na.rm=T)

write(paste('Figure S4 results: BvPrP RT-QuIC sensitivity to positive controls was: ',percent(poscon_positives/(poscon_positives+poscon_negatives)),'\n',sep=''),text_stats_path,append=T)


