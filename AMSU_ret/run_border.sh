# important directories:
BASE_PATH = /freax/storage/home/pmills/isoline_ret2
VEC_PATH = ${BASE_PATH}/data/training
BIN_PATH = $(BASE_PATH)/bin

# zenith angle:
ZANGLE = 0.00

# surface type:
SURFACE = land

# retrieval parameters:
SIGMA = 36
SQ = 0.001

# classification parameters:
NNN = 1000
WC = 30
NSAMPLES = 1000
TOL = 0.0001
VAR0 = 4

VECFILE = $(VEC_PATH)/est_bt_cld_$(SURFACE).vec

# base name for class borders:
BORDER_BASE=iso$(ZENITH_ANGLE)z$(SIGMA)s$(SQ)sq_$(SURFACE)

#file containing wv vmrs for training data:
TRAIN_SQ = train_SQ$(SQ).dat

# here we make the test dates:
YEAR = 2003
MONTH = 1
DAYLIST = 1 2 3 4 5
TEST_DAYS = $(foreach DAY, $(DAYLIST), $(YEAR)/$(MONTH)/$(DAY))
HOURS = 0 6 12 18
TEST_DATES = $(foreach HOUR, $(HOURS), $(TEST_DAYS)-$(HOUR))

${BIN_PATH}/get_date_inds ${PROF_FILE} ${TEST_DATES} > test_inds.txt
${BIN_PATH}/get_sigma_level3 $(PROF_FILE) train_SQ $(SIGMA) > dum.txt
${BIN_PATH}/prune_vecfile ${VECFILE} train.vec < test_inds.txt
${BIN_PATH}/float_to_class ${TRAIN_SQ} ${SQ} train1.cls
${BIN_PATH}/prune_clsfile train1.cls train.cls < test_inds.txt
time multi_borders -k ${NNN} -w ${WC} -s ${NSAMPLES} -t ${TOL} \
			-v ${VAR0} train ${BORDER_BASE}

# for pruning the test dates from the training data:
test_inds.txt: $(PROF_FILE) 

$(TRAIN_SQ): $(PROF_FILE)


