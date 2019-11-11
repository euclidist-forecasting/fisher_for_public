# Instructions for using the code in this repo

Generate a header for the top of your file for example:

```bash
python genhead.py -p wm wb h ns -d lnH lnDa lnfs8 lnbs8 Ps
```
Use the `-f` flag to input a file to be re-written!

Make a Fisher matrix density plot for the shape parameters only for e.g.:

```bash
python plfmat.py [inputfiles] -c 5 -p wm wb h ns -f
```

Make marginalised and unmarginalised error plots for shape parameters only for e.g.:

```bash
python plfmat.py [inputfiles] -c 5 -p wm wb h ns
```

`fishlib.py` and `densplotfisher.py` are only libraries for some of the functions used in the plotting.
