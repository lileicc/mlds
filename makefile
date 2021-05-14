MLDSFILE=MLE_mlds.m \
  initialize_parameters.m \
  number_of_parameters.m \
  traceprod.m \
  learn_mlds.m \
  readme.txt \
  update_multilinear_operator.m \
  demo_synthetic.m \
  demo_sst.m \
  license.txt \
  set_optional_argument.m \
  vec.m \
  descend.m \
  makefile \
  subcell.m \
  vec2ten.m \
  mat2ten.m \
  ten2mat.m \
  err.m \
  mkron.m \
  ten2vec.m

DYNAMMOFILE=dynammo

DATAFILE=data

REVNUM:=$(shell svn info |grep Revision: |cut -c11-)

VERSION:=$(shell echo "scale=1;(${REVNUM} - ${REVNUM} % 100) / 1000" | bc)

demo:
	matlab -r demo_synthetic

zip: tar

tar: ${DYNAMMOFILE} ${MLDSFILE} 
	mkdir -p temp/mlds
	cp -r ${DYNAMMOFILE} temp/mlds/
	cp ${MLDSFILE} temp/mlds/
	cp -r ${DATAFILE} temp/mlds/
	cd temp; zip -r mlds-r${REVNUM}.zip mlds
	cp temp/mlds-r${REVNUM}.zip ./
	rm -r -f temp

sync: 
	svn up
	svn commit -m "commit" 
	svn up
