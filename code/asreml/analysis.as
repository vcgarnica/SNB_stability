####################################################################
######################## Vinicius GxE Study ########################

!ARGS 8 t50 !RENAME 2 !WORKSPACE 20 !OUTFOLDER C:\Users\mshaliz\Desktop\Vinicius\analysis3\out

Title: Stability of wheat to disease

 site  !A
 block  !A
 cult  !A
 plot  !A
 trt  !A
 eu  !A
 sev
 audps
 raudps
 naudps
 omega !*100
 t50


!FOLDER C:\Users\mshaliz\Desktop\Vinicius\analysis3
# data
dat.csv !SKIP 1 !DOPART $A !MVINCLUDE !MP 15

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# single site analysis
# trait = sev
!PART 1
!CONTINUE !MSV !NODISPLAY
$B ~ mu at(site,CL22,RW22,SB22,MR22,PY22,KS22,OX23,KS23,PY23,UN23,SB23,LB23,AL24,BE24,KS24,LB24,RO24,SB24).block,
     !r diag(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# single site analysis
# trait = raudps
!PART 2
!CONTINUE !MSV !NODISPLAY
$B ~ mu at(site,CL22,RW22,SB22,MR22,PY22,KS22,OX23,KS23,PY23,UN23,SB23,LB23,AL24,BE24,KS24,LB24,RO24,SB24).block,
     !r diag(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# single site analysis
# trait = omega
!PART 3
!CONTINUE !MSV !NODISPLAY
$B ~ mu at(site,CL22,RW22,SB22,MR22,PY22,KS22,OX23,KS23,PY23,UN23,SB23,LB23,AL24,BE24,KS24,LB24,RO24,SB24).block,
     !r diag(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# single site analysis
# trait = t50
!PART 4
!CONTINUE !MSV !NODISPLAY
$B ~ mu at(site,CL22,RW22,SB22,MR22,PY22,KS22,OX23,KS23,PY23,UN23,SB23,LB23,AL24,BE24,KS24,LB24,RO24,SB24).block,
     !r diag(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# multi site analysis
# trait = sev
!PART 5
!CONTINUE !MSV !NODISPLAY
$B ~ mu site site.block,
     !r xfa3(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE

predict cult
predict site cult


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# multi site analysis
# trait = raudps
!PART 6
!CONTINUE !MSV !NODISPLAY
$B ~ mu site site.block,
     !r xfa3(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE

predict cult
predict site cult


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# multi site analysis
# trait = omega
!PART 7
!CONTINUE !MSV !NODISPLAY
$B ~ mu site site.block,
     !r xfa2(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE

predict cult
predict site cult


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# multi site analysis
# trait = t50
!PART 8
!CONTINUE !MSV !NODISPLAY
$B ~ mu site site.block,
     !r xfa2(site).cult
        residual sat(site).id(units)
VPREDICT !DEFINE

predict cult
predict site cult









