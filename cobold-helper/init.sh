#!/bin/bash
# Set up a new run provided an ID
# Just copies all the files from outer directories to 'runs'
# Assumes the start, par files exist for the given ID
# Copies 'template.{cmd,job}' and renames template > model ID
# Does NOT change 'chem.dat', this must still
# be edited manually for chem runs
if [ $# -eq 0 ]
then
	echo "Error: Must provide a model ID, e.g. d2tg25mm00"
	exit
fi

COBOLD_DIR=/home/sd/cobold
RUN_DIR=$COBOLD_DIR/runs
model_id=$1

# Check for 'exe' option
if [ $# -eq 1 ]; then  # no 'exe' option provided
	exe=$RUN_DIR/rhd.exe
fi

# Should check exe via flag, not position
if [ $# -ge 2 ]; then
	exe=$COBOLD_DIR/exe/$2
fi

# Copy files and rename templates
# Start model and parameter file must exist
cp $RUN_DIR/start/$model_id.start ./rhd.start
cp $COBOLD_DIR/par/$model_id.par $RUN_DIR/rhd.par

# Copy templates
cp $COBOLD_DIR/cmd/template.cmd $RUN_DIR/cobold.cmd
cp $COBOLD_DIR/job/template.job $RUN_DIR/cobold.job

# Edit template cmd and job file, replacing 'template' with 'model_id'
old="template"
sed -i "s|$old|$model_id|g;" $RUN_DIR/cobold.cmd
sed -i "s|$old|$model_id|g;" $RUN_DIR/cobold.job

# Exe if specified
cp $exe $RUN_DIR/rhd.exe

# Check for --single option to edit 'cmd' file for 1 run
while [ ! $# -eq 0 ]
do
	case "$1" in
		--single | -s)
			sed -ni 's/doc/eoc/;1p;' $RUN_DIR/cobold.cmd
			;;
	esac
	shift
done

# Create directories if necessary
mkdir -p {$RUN_DIR/full,$RUN_DIR/mean,$RUN_DIR/fine,$RUN_DIR/out}/$model_id

# Print to terminal to make sure it worked
find . -name $model_id -type d 
head -1 cobold.cmd
tail -2 cobold.cmd
grep 'IDDIR' cobold.job | head -1

echo "Check above statements to see if it's correct!"
echo "If using CHEM, ensure the supplied start model has correct abundances!"
echo "If using CHEM, ensure 'chem.dat' is the correct network!"
echo "Provide correct 'rhd.exe' from 'exe' dir if you want to change it!"
