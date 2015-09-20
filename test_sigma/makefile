# data paths:
BASE_PATH = /home/petey/my_software/isoline_ret2
VECFILE_PATH = $(BASE_PATH)/data/testfields/vecfiles/Wyser
TRAINVEC_PATH = $(BASE_PATH)/data/training
FIELDFILE_PATH = $(BASE_PATH)/data/testfields
WORK_PATH = $(BASE_PATH)/data/work

# system commands:
GV = gsview32

# commands:
RTSIM_TRAIN = /freax/storage/home/pmills/rtsim_bremen2/rtsim_train8

# retrieval paramters:
#DATE = 1998082718
DATE = 2000101600
#DATE = 2003010100
SQ = 0.001
SIGMA = 36

ZENITH_ANGLE=0.00

# classification parameters:
TRAINING=.b
#NNN=200
NNN=-1
WC=30
NSAMPLES=250
TOL=0.0001
VAR0=4
NORM_FLAG = -n
#glue the parameters together (most of them):
KERNEL_PARM = - k $(NNN) -w $(WC) -v $(VAR0)

# cloudy or clearsky:
CLD=cld

# plotting parameters:

# filenames:
# base-names of the classification retrievals:
FBASE = ecmwf$(DATE)_$(CLD)$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z
FBASE_CLD = ecmwf$(DATE)_cld$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z

# files containing simulated bt vectors for ecmwf fields:
VECFILE_CLD = $(VECFILE_PATH)/ecmwf$(DATE)_cld$(ZENITH_ANGLE)z.vec
VECFILE_CLD_SEA = $(VECFILE_PATH)/ecmwf$(DATE)_cld_sea$(ZENITH_ANGLE)z.vec

# files containing training data:
TRAIN_VEC = $(TRAINVEC_PATH)/est_bt_cld_land_z$(ZENITH_ANGLE).vec
TRAIN_VEC_SEA = $(TRAINVEC_PATH)/est_bt_cld_sea_z$(ZENITH_ANGLE).vec
PROFILE_FILE = ../data/est_ecmwf_profiles.dat
TRAIN_SQ = train_SQ$(SIGMA)s.dat

# all possible dates, in the more usual format:
ALL_DATES = 1998/8/27-18 1999/3/2-18 1999/9/22-12 2000/10/16-0 2001/12/20-6 2002/6/12-12

SURFACE =

# for storing brightness temperature calculations:
BT_PATH = $(BASE_PATH)/data/testfields/bts/Wyser

VPATH = $(BT_PATH)

all: $(FBASE).ret$(TRAINING).cls $(FBASE).ret$(TRAINING).con

all_dates:
	make DATE=1998082718
	make DATE=1999030218
	make DATE=1999092212
	make DATE=2000101600
	make DATE=2001122006
	make DATE=2002061212

help:
	@echo
	@echo "SUMMARY:"
	@echo
	@echo "make -f makefile.ret [all] [check] [comp] [plot] [help]"
	@echo "     [clean_ret] [clean_train] [clean_all] " 
	@echo "     [DATE = YYYYMMDDHH] [PRESSURE = \<plev\>] [SQ = <SQ>]"
	@echo "     [TRAINING = [.b]] [NNN = <k>] [WC = <Wc>]" 
	@echo "     [NSAMPLES = <nsamples>] [TOL = <tol>]"
	@echo "     [CLD = [cld]] [SURFACE = [land | sea | ]]"
	@echo
	@echo "This makefile is designed to test the isoline retrieval algorithm."
	@echo "It accepts the following targets:"
	@echo
	@echo "all:     Performs the retrieval for default or specified parameters."
	@echo "check:   Checks the accuracy of the retrieval."
	@echo "comp:    Compares retrievals without and with training."
	@echo "plot:    Plots the retrieval."
	@echo "help:    Displays this help screen."
	@echo
	@echo "clean_ret:   Cleans out the classification retrievals."
	@echo "clean_train: Cleans out the pre-trained border files."
	@echo 
	@echo "...and the following macros:"
	@echo 
	@echo "Retrieval parameters:"
	@echo "DATE      = Date of desired retrieval in the following format: YYYYMMDDHH"
	@echo "SIGMA     = Sigma level on which to perform retrieval."
	@echo "SQ        = Mass-mixing-ratio (specific humidity) of desired isoline."
	@echo 
	@echo "Classifier algorithm parameters:"
	@echo "TRAINING  = Set to \"b.\" to specify a separate training phase."
	@echo "            This is the default.  To leave out training phase, set to null."
	@echo "NNN       = Number of nearest neighbours to use in classification scheme."
	@echo "WC        = Objective critical weight."
	@echo "NSAMPLES  = Number of times to sample the border (if there is a training phase.)"
	@echo "TOL       = Tolerance of each border sample."
	@echo
	@echo "*****  Note:  changes in any but the first classifier parameter do NOT result"
	@echo "       in updates to the relevant files.  Files must be deleted by hand."
	@echo 
	@echo "Simulation parameters:"
	@echo "CLD       = Set to \"cld\" to include clouds and scattering in both"
	@echo "            the test and training data.  This is the default."
	@echo "            To perform clearsky retrievals, leave blank.  Mixed retrievals"
	@echo "            are not yet possible."
	@echo "SURFACE   = Surface type for retrieval.  \"land\" assumes all surface is land,"
	@echo "            \"sea\" assumes only sea surface.  If set to null, performs a"
	@echo "            mixed retrieval.  This is the default."
	@echo
	@echo "These are the current defaults:"
	@echo "DATE     = $(DATE)"
	@echo "SIGMA    = $(SIGMA)"
	@echo "SQ       = $(SQ)"
	@echo "NNN      = $(NNN)"
	@echo "WC       = $(WC)"
	@echo "NSAMPLES = $(NSAMPLES)"
	@echo "TOL      = $(TOL)"

clean_all: clean_ret clean_train
	rm -f *.ps
	rm -f *.dat

# clean out retrievals:
clean_ret:
	rm -f *.cls
	rm -f *.con

clean_plot:
	rm -f *.ps

# clean out border files:
clean_train:
	rm -f *.brd
	rm -f *.bgd
	rm -f *.std

# check the accuracy:
check: ecmwf$(DATE)_$(SIGMA)s$(SQ)SQ.cls \
		$(FBASE).$(SURFACE)ret$(TRAINING).con \
		$(FBASE).$(SURFACE)ret$(TRAINING).cls
	cls_comp_stats ecmwf$(DATE)_$(SIGMA)s$(SQ)SQ.cls \
			$(FBASE).$(SURFACE)ret$(TRAINING)

# compare with training and without:
comp: $(FBASE).$(SURFACE)ret.cls $(FBASE).$(SURFACE)ret.b.cls
	cls_comp_stats $(FBASE).$(SURFACE)ret.cls $(FBASE).$(SURFACE)ret.b.cls

# plot the postscript file:
plot: $(FBASE).$(SURFACE)plot$(TRAINING).ps
	$(GV) $(FBASE).$(SURFACE)plot$(TRAINING).ps &

# this creates the postscript file for plotting:

ORIENT = -ZTLa

$(FBASE).$(SURFACE)plot$(TRAINING).ps: $(FBASE).$(SURFACE)ret$(TRAINING).cls \
		$(FBASE).$(SURFACE)ret$(TRAINING).con \
		ecmwf$(DATE)SQ$(SIGMA)s.dat makefile
	echo $(SQ) A > SQ_thresh.txt
	../bin/print_bin_file $(FBASE).ret$(TRAINING).con | xyz2grd -R0/358.5/-90/90 \
		-I1.5/1.5 $(ORIENT) -Gcon.grd
	../bin/print_bin_file ecmwf${DATE}SQ$(SIGMA)s.dat | \
		xyz2grd -R0/358.5/-90/90 -I1.5/1.5 $(ORIENT) -GSQ.grd
	psbasemap -R0/360/-90/90 -JQ180/18 -Bf15a30g30 -K > $@
	grdimage con.grd -R -J -Cconfidence.cpt -O -K >> $@
	pscoast -R -J -Dl -W -O -K >> $@
	#grdcontour SQ.grd -R -J -CSQ_aux.txt -A0.003 -O -K \
	#	 -Wa0.5p,black,solid >> $@
	grdcontour SQ.grd -R -J -CSQ_thresh.txt -O -K \
		-Wa1.5p,black,solid >> $@
	../bin/ecmwf_cls_field_to_xyz $(FBASE).ret$(TRAINING).cls | psxy -R -J -O -Sc0.01 >> $@

# this creates the true class values for the given date, pressure level and SQ theshold:
ecmwf$(DATE)_$(SIGMA)s$(SQ)SQ.cls: ecmwf$(DATE)SQ$(SIGMA)s.dat
	../bin/float_to_class  ecmwf$(DATE)SQ$(SIGMA)s.dat $(SQ) $@

# this extracts the volume mixing ratios for the given pressure level:
ecmwf$(DATE)SQ$(SIGMA)s.dat: $(FIELDFILE_PATH)/ecmwf$(DATE).dat
	../bin/get_sigma_level3 -s $(FIELDFILE_PATH)/ecmwf$(DATE).dat ecmwf$(DATE)SQ $(SIGMA) > dum.txt

# here we unite the retrieval over land and retrieval over sea:
$(FBASE).ret$(TRAINING).cls: $(FBASE).landret$(TRAINING).cls \
		$(FBASE).searet$(TRAINING).cls landmask.txt
	../bin/unite_landsea_cls landmask.txt $(FBASE).landret$(TRAINING).cls \
			$(FBASE).searet$(TRAINING).cls $@

# here we unite the retrieval over land and retrieval over sea:
$(FBASE).ret$(TRAINING).con: $(FBASE).landret$(TRAINING).con \
		$(FBASE).searet$(TRAINING).con landmask.txt
	../bin/unite_landsea_cls landmask.txt $(FBASE).landret$(TRAINING).con \
			$(FBASE).searet$(TRAINING).con $@

# this extracts the volume-mixing-ratios for the training data:
train_SQ$(SIGMA)s.dat:  $(PROFILE_FILE)
	../bin/get_sigma_level3 $(PROFILE_FILE) train_SQ $(SIGMA) > dum.txt


################     cloudy retrievals     #################################

# do the retrieval over land:
$(FBASE_CLD).landret.b.cls $(FBASE_CLD).landret.b.con: $(VECFILE_CLD) \
		cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.brd \
		cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.bgd
	classify_b $(NORM_FLAG) -n cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z \
			$(VECFILE_CLD) $(FBASE_CLD).landret.b

# do the retrieval over sea:
$(FBASE_CLD).searet.b.cls $(FBASE_CLD).searet.b.con: $(VECFILE_CLD_SEA) \
		cld_sea_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.brd \
		cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.brd
	classify_b $(NORM_FLAG) cld_sea_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z \
			$(VECFILE_CLD_SEA) $(FBASE_CLD).searet.b

# base name for class borders:
#BORDER_BASE_LAND=../ret/iso$(ZENITH_ANGLE)z$(SIGMA)s$(SQ)sq_land
#BORDER_BASE_SEA=../ret/iso$(ZENITH_ANGLE)z$(SIGMA)s$(SQ)sq_sea

# do the retrieval over land:
#$(FBASE_CLD).landret.b.cls $(FBASE_CLD).landret.b.con: $(VECFILE_CLD) \
#		$(BORDER_BASE_LAND).brd $(BORDER_BASE_LAND).bgd
#	classify_b $(BORDER_BASE_LAND) $(VECFILE_CLD) $(FBASE_CLD).landret.b

# do the retrieval over sea:
#$(FBASE_CLD).searet.b.cls $(FBASE_CLD).searet.b.con: $(VECFILE_CLD_SEA) \
#		$(BORDER_BASE_SEA).brd $(BORDER_BASE_SEA).bgd
#	classify_b $(BORDER_BASE_SEA) $(VECFILE_CLD_SEA) $(FBASE_CLD).searet.b

# generate the class borders for the retrieval over land:
cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.brd \
		cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.bgd: \
		$(TRAIN_VEC) $(TRAIN_SQ) test_inds.txt
	../bin/prune_vecfile $(TRAIN_VEC) train.vec < test_inds.txt
	../bin/float_to_class $(TRAIN_SQ) $(SQ) train1.cls
	../bin/prune_clsfile train1.cls train.cls < test_inds.txt
	time class_borders $(NORM_FLAG) $(KERNEL_PARM) -s $(NSAMPLES) -t $(TOL) -v \
			train cld_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z

# generate the class borders for the retrieval over sea:
cld_sea_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.brd \
		cld_sea_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z.bgd: \
		$(TRAIN_VEC_SEA) $(TRAIN_SQ) test_inds.txt
	../bin/prune_vecfile $(TRAIN_VEC_SEA) train_sea.vec < test_inds.txt
	../bin/float_to_class $(TRAIN_SQ) $(SQ) train1_sea.cls
	../bin/prune_clsfile train1_sea.cls train_sea.cls < test_inds.txt
	time class_borders $(NORM_FLAG) $(KERNEL_PARM) -s $(NSAMPLES) -t $(TOL) -v \
			train_sea cld_sea_iso$(SIGMA)s$(SQ)SQ$(ZENITH_ANGLE)z


#################     Direct retrievals      ###########################

################     cloudy:           #################################
# do the retrieval over land:
$(FBASE_CLD).landret.cls $(FBASE_CLD).landret.con: $(VECFILE_CLD) \
		$(TRAIN_SQ) $(TRAIN_VEC_SEA)
	../bin/prune_vecfile $(TRAIN_VEC_SEA) train_sea.vec < test_inds.txt
	../bin/float_to_class $(TRAIN_SQ) $(SQ) train1.cls
	../bin/prune_clsfile train1.cls train.cls < test_inds.txt
	time classify_a $(NORM_FLAG) $(KERNEL_PARM) train $(VECFILE_CLD) $(FBASE_CLD).landret 
	rm train.vec train.cls

# do the retrieval over sea:
$(FBASE_CLD).searet.cls $(FBASE_CLD).searet.con: $(VECFILE_CLD_SEA) \
		$(TRAIN_SQ) $(TRAIN_VEC)
	../bin/prune_vecfile $(TRAIN_VEC) train.vec < test_inds.txt
	../bin/float_to_class $(TRAIN_SQ) $(SQ) train1_sea.cls
	../bin/prune_clsfile train1_sea.cls train_sea.cls < test_inds.txt
	time classify_a $(NORM_FLAG) $(KERNEL_PARM) train_sea $(VECFILE_CLD_SEA) $(FBASE_CLD).searet
	rm train_sea.vec train_sea.cls


#############   Removing the test dates from the training data:    ########

#train.vec: $(TRAIN_VEC) test_inds.txt
#	../bin/prune_vecfile $(TRAIN_VEC) train.vec < test_inds.txt

#train_sea.vec: $(TRAIN_VEC_SEA) test_inds.txt
#	../bin/prune_vecfile $(TRAIN_VEC_SEA) train_sea.vec < test_inds.txt

test_inds.txt: $(PROFILE_FILE)
	../bin/get_date_inds $(PROFILE_FILE) $(ALL_DATES) > test_inds.txt

#############    For making the .vec files themselves:             ########

#$(VECFILE_CLD): bt$(DATE)A$(ZENITH_ANGLE)z.dat bt$(DATE)B$(ZENITH_ANGLE)z.dat \
#		../bin/rttov_sim_collect
#	../bin/rttov_sim_collect $(BT_PATH)/bt$(DATE)A$(ZENITH_ANGLE)z.dat \
#			$(BT_PATH)/bt$(DATE)B$(ZENITH_ANGLE)z.dat \
#			$(VECFILE_CLD)

bt$(DATE)A$(ZENITH_ANGLE)z.dat: $(FIELDFILE_PATH)/ecmwf$(DATE).dat
	cd $(WORK_PATH)
	$(RTSIM_TRAIN) $(FIELDFILE_PATH)/ecmwf$(DATE).dat \
			bt$(DATE)A$(ZENITH_ANGLE)z.dat 16 3 1 $(ZENITH_ANGLE) 1000
	mv bt$(DATE)A$(ZENITH_ANGLE)z.dat $(BT_PATH)

bt$(DATE)B$(ZENITH_ANGLE)z.dat: $(FIELDFILE_PATH)/ecmwf$(DATE).dat
	cd $(WORK_PATH)
	$(RTSIM_TRAIN) $(FIELDFILE_PATH)/ecmwf$(DATE).dat \
			bt$(DATE)B$(ZENITH_ANGLE)z.dat 16 4 1 $(ZENITH_ANGLE) 1000
	mv bt$(DATE)B$(ZENITH_ANGLE)z.dat $(BT_PATH)

#$(VECFILE_CLD_SEA): bt$(DATE)A_sea$(ZENITH_ANGLE)z.dat \
#		bt$(DATE)B_sea$(ZENITH_ANGLE)z.dat ../bin/rttov_sim_collect
#	../bin/rttov_sim_collect $(BT_PATH)/bt$(DATE)A_sea$(ZENITH_ANGLE)z.dat \
#			$(BT_PATH)/bt$(DATE)B_sea$(ZENITH_ANGLE)z.dat \
#			$(VECFILE_CLD_SEA)

bt$(DATE)A_sea$(ZENITH_ANGLE)z.dat: $(FIELDFILE_PATH)/ecmwf$(DATE).dat
	cd $(WORK_PATH)
	$(RTSIM_TRAIN) $(FIELDFILE_PATH)/ecmwf$(DATE).dat \
			bt$(DATE)A_sea$(ZENITH_ANGLE)z.dat 16 3 0 $(ZENITH_ANGLE) 1000
	mv bt$(DATE)A_sea$(ZENITH_ANGLE)z.dat $(BT_PATH)

bt$(DATE)B_sea$(ZENITH_ANGLE)z.dat: $(FIELDFILE_PATH)/ecmwf$(DATE).dat
	cd $(WORK_PATH)
	$(RTSIM_TRAIN) $(FIELDFILE_PATH)/ecmwf$(DATE).dat \
			bt$(DATE)B_sea$(ZENITH_ANGLE)z.dat 16 4 0 $(ZENITH_ANGLE) 1000
	mv bt$(DATE)B_sea$(ZENITH_ANGLE)z.dat $(BT_PATH)

$(FIELDFILE_PATH)/ecmwf$(DATE).dat: $(FIELDFILE_PATH)/ecmwf$(DATE)_proflist.txt
	../bin/make_ecmwf_profile_file $(FIELDFILE_PATH)/ecmwf$(DATE)_proflist.txt \
			$(FIELDFILE_PATH)/ecmwf$(DATE).dat

$(FIELDFILE_PATH)/ecmwf$(DATE)_proflist.txt:
	../bin/ecmwf_field_indices $(DATE) > $(FIELDFILE_PATH)/ecmwf$(DATE)_proflist.txt

