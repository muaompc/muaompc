VERSION=1.0.0-b0
prb=mpc_myprb
find .. | grep -E "(__pycache__|\.pyc|\.pyo$)" | xargs rm -rf
cd ../examples/ldt
cd arduino
rm -rf muaompc_arduino ${prb} ${prb}_mydat
cd ..
cd matlab_interface
rm -rf ${prb}
cd ..
cd path_following_cython
rm -rf ${prb}
cd ..
cd trajtrack
rm -rf ${prb}
cd ..
cd tutorial
rm -rf ${prb}
cd ..
cd tutorial_advanced
rm -rf ${prb}
cd ../../..
rm -rf mpc_test*
rm -rf fixtures/mpc_test
echo 'Clean up ready'
rm -rf build
rm -rf dist
python -m build && pip install ~/github/muaompc/muaompc/dist/muaompc-${VERSION}-py3-none-any.whl --force
echo 'Distribuition ready'
echo 'muaompc install ready to generate new documentation'
cd doc
make clean
make latexpdf
make html
cd ..
echo 'Generation of Docu ready'

cd doc/_build/
tar czf html.tar.gz html
cp html.tar.gz  ../muaompc-${VERSION}.tar.gz
cp latex/muaompc.pdf  ../muaompc-${VERSION}.pdf
