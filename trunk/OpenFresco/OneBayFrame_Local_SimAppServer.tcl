# File: OneBayFrame_Local_SimAppServer.tcl (use with OneBayFrame_Local_Client.tcl)
#
# $Revision: $
# $Date: $
# $URL: $
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
# Created: 11/06
# Revision: A
#
# Purpose: this file contains the tcl input to perform
# a local hybrid simulation of a one bay frame
# with two experimental twoNodeLink elements.
# The specimen is simulated using the SimUniaxialMaterials
# controller.


# ------------------------------
# Start of model generation
# ------------------------------
# create ModelBuilder (with two-dimensions and 2 DOF/node)
model BasicBuilder -ndm 2 -ndf 2

# Define materials
# ----------------
# uniaxialMaterial ElasticForce $matTag $E
uniaxialMaterial ElasticForce 1 2.8
#uniaxialMaterial Elastic 1 2.8

# Define experimental control
# ---------------------------
# expControl SimUniaxialMaterialsForce $tag $matTags
#expControl SimUniaxialMaterials 1 1
expControl SimUniaxialMaterialsForce 1 1
# expControl xPCtargetForce $tag ipAddr $ipPort 
#expControl xPCtargetForce 1 "192.168.2.20" 22222 HybridControllerF1F1_Filter_1Act "D:/PredictorCorrector/RTActualTestModels/cmAPI-xPCTarget-STS"
#expControl xPCtarget 1 1 "192.168.2.20" 22222 HybridControllerD1D1_Filter_1Act "D:/PredictorCorrector/RTActualTestModels/cmAPI-xPCTarget-STS"

# Define experimental setup
# -------------------------
# expSetup OneActuator $tag <-control $ctrlTag> $dir -sizeTrialOut $t $o <-trialDispFact $f> ...
expSetup OneActuator 1 -control 1 1 -sizeTrialOut 1 1

# Define experimental site
# ------------------------
# expSite LocalSite $tag $setupTag
expSite LocalSite 1 1

# ------------------------------
# End of model generation
# ------------------------------


# ------------------------------
# Start of recorder generation
# ------------------------------
# create the recorder objects
expRecorder Control -file Control_ctrlFrc.out -time -control 1 ctrlForce
# --------------------------------
# End of recorder generation
# --------------------------------


# ------------------------------
# Start the server process
# ------------------------------
# startSimAppSiteServer $siteTag $port <-ssl>
startSimAppSiteServer 1 8090;  # use with experimental element in FEA
# --------------------------------
# End of analysis
# --------------------------------
