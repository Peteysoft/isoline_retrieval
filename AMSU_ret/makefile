# important directories:
BASE_PATH = /freax/storage/home/pmills/isoline_ret2
VEC_PATH = $(BASE_PATH)/data/training
BIN_PATH = $(BASE_PATH)/bin

AMSU1C_PATH = $(BASE_PATH)/data/amsul1c

# retrieval parameters:
SIGMA = 36
SQ = 0.001

# classification parameters:
NNN = 3000
WC = 30
NSAMPLES = 500
TOL = 0.0001
VAR0 = 4

# file listing all profiles in the training data:
PROF_LIST = $(BASE_PATH)/data/est_prof_list

# data file containing all the profiles:
PROF_FILE = $(BASE_PATH)/data/est_ecmwf_profiles.dat

# file containing vector training data:
VECFILE = $(VEC_PATH)/est_bt_cld_$(SURFACE)_z$(ZENITH_ANGLE).vec

# base name for class borders:
BORDER_BASE=iso$(ZENITH_ANGLE)z$(SIGMA)s$(SQ)sq_$(SURFACE)

# initialization file:
INITFILE = isoret.init

#file containing wv vmrs for training data:
TRAIN_SQ = train_SQ$(SIGMA)s.dat

# all the zenith angles:
ZANGLES = 0.00 2.80 6.17 9.55 12.92 16.29 19.66 23.03 26.40 29.77 33.14 36.52 \
		39.89 43.26 46.63 50.00 60.02

# here we make the test dates:
YEAR = 2002
MONTH = 9
DAYLIST = 1 2 3 4 5 6 7 8 9
TEST_DAYS = $(foreach DAY, $(DAYLIST), $(YEAR)/$(MONTH)/$(DAY))
HOURS = 0 6 12 18
TEST_DATES = $(foreach HOUR, $(HOURS), $(TEST_DAYS)-$(HOUR))

# list of files containing class borders for each zenith angle:
BORDER_FILES = $(foreach SURFACE, land sea, $(foreach ZENITH_ANGLE, $(ZANGLES), $(BORDER_BASE).0000.bgd))

# list of amsu 1c files:
AMSUB_DIRS = $(foreach DAY, $(DAYLIST), $(AMSU1C_PATH)/amsub$(YEAR)$(MONTH)$(DAY))
AMSUB_FILES=$(foreach DIR, $(AMSUB_DIRS), $(DIR)/$(shell ls $(DIR)))

AMSUA1C = $(patsubst AMBX, AMAX, $(AMSUB1C))

RETFILE = $(patsubst NSS\.AMBX\.N?, RET$(SIGMA)s$(SQ)sq, $(base $(AMSUB1C))).r1c

$(RETFILE): $(AMSUB1C) $(AMSUA1C) $(BORDER_FILES) $(INITFILE)
	$(BIN_PATH)/retrieve_isoline $(AMSUA1C) $(AMSUB1C) $(RETFILE)

$(INITFILE):
	echo 6 > $(INITFILE)
	echo -e "$(foreach ZENITH_ANGLE, $(ZANGLES), $(ZENITH_ANGLE)\n $(foreach SURFACE, land sea, $(BORDER_BASE)\n))" >> $(INITFILE)

border_files:
	make ZENITH_ANGLE=0.00 SURFACE=sea border
	make ZENITH_ANGLE=0.00 SURFACE=land border
	make ZENITH_ANGLE=2.80 SURFACE=sea border
	make ZENITH_ANGLE=2.80 SURFACE=land border
	make ZENITH_ANGLE=6.17 SURFACE=sea border
	make ZENITH_ANGLE=6.17 SURFACE=land border
	make ZENITH_ANGLE=9.55 SURFACE=sea border
	make ZENITH_ANGLE=9.55 SURFACE=land border
	make ZENITH_ANGLE=12.92 SURFACE=sea border
	make ZENITH_ANGLE=12.92 SURFACE=land border
	make ZENITH_ANGLE=16.29 SURFACE=sea border
	make ZENITH_ANGLE=16.29 SURFACE=land border
	make ZENITH_ANGLE=19.66 SURFACE=sea border
	make ZENITH_ANGLE=19.66 SURFACE=land border
	make ZENITH_ANGLE=23.03 SURFACE=sea border
	make ZENITH_ANGLE=23.03 SURFACE=land border
	make ZENITH_ANGLE=26.40 SURFACE=sea border
	make ZENITH_ANGLE=26.40 SURFACE=land border
	make ZENITH_ANGLE=29.77 SURFACE=sea border
	make ZENITH_ANGLE=29.77 SURFACE=land border
	make ZENITH_ANGLE=33.14 SURFACE=sea border
	make ZENITH_ANGLE=33.14 SURFACE=land border
	make ZENITH_ANGLE=36.52 SURFACE=sea border
	make ZENITH_ANGLE=36.52 SURFACE=land border
	make ZENITH_ANGLE=39.89 SURFACE=sea border
	make ZENITH_ANGLE=39.89 SURFACE=land border
	make ZENITH_ANGLE=43.26 SURFACE=sea border
	make ZENITH_ANGLE=43.26 SURFACE=land border
	make ZENITH_ANGLE=46.63 SURFACE=sea border
	make ZENITH_ANGLE=46.63 SURFACE=land border
	make ZENITH_ANGLE=50.00 SURFACE=sea border
	make ZENITH_ANGLE=50.00 SURFACE=land border
	make ZENITH_ANGLE=60.02 SURFACE=sea border
	make ZENITH_ANGLE=60.02 SURFACE=land border

border: $(BORDER_BASE).brd $(BORDER_BASE).brd $(BORDER_BASE).std

$(BORDER_FILES):
	$(foreach SURFACE, land sea, $(foreach ZENITH_ANGLE, $(ZANGLES), \
	make ZENITH_ANGLE=$(ZENITH_ANGLE) SURFACE=$(SURFACE) border))

$(BORDER_BASE).brd $(BORDER_BASE).bgd \
		$(BORDER_BASE).std: $(VECFILE) $(TRAIN_SQ) test_inds.txt
	$(BIN_PATH)/prune_vecfile $(VECFILE) train.vec < test_inds.txt
	$(BIN_PATH)/float_to_class $(TRAIN_SQ) $(SQ) train1.cls
	$(BIN_PATH)/prune_clsfile train1.cls train.cls < test_inds.txt
	time class_borders -n -k $(NNN) -w $(WC) -s $(NSAMPLES) -t $(TOL) \
			-v $(VAR0) train $(BORDER_BASE)

# for pruning the test dates from the training data:
test_inds.txt: $(PROF_FILE) makefile
	$(BIN_PATH)/get_date_inds $(PROF_FILE) $(TEST_DATES) > test_inds.txt

$(TRAIN_SQ): $(PROF_FILE)
	$(BIN_PATH)/get_sigma_level3 $(PROF_FILE) train_SQ $(SIGMA) > dum.txt


