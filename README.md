Nathan Gillispie's knot theory repository.

Due to the large number of PDF and Mathematica Notebook files, this is only a repository for backing up my work rather than an actively maintained code base.

The projects of interest here are Relax_3D, Knot Relaxer, PDToConway.nb, and NewDrawPD-explained.nb

## Relax_3D

Relax_3D is my first ever knot theory related project. It is written in Processing (Java), at a time when I was trying to understand knot energy functions. I was fascinated with the concept, and thought it would be interesting to make software that minimized this energy. Unfortunately, I had only taken calculus 1 and many of the functions operated on vector calculus, which I had to develop an intuition for, myself. It was very exciting to get working in the 2D case and finally here, the 3D case.

## Knot_Relaxer

CPS is an infamous course (love it or hate it) that all Gatton students take their second semester. I decided, since in knot theory we typically work in Mathematica, that the high level code would be easier to read and add menus. So Relax_3D was reimplemented here.

## PD_Conway

An ongoing concern with knot theory is the many many ways there are to represent knots in text form. Some of them are good for one thing but not others. For example, Planar Diagram (PD) code and Dowker-Thistlewaite (DT) code are best for graphical representation, however, Conway notation is not. When reading conway notation, it can be difficult to visualize the results, so ConwayToPD was made for this reason. PD code however, is what our research group used for all our knot simplification algorithms, therefore, to convert back to Conway notation, PDtoConway was made. This is not a trivial task however, as all conway polyhedra have to be hard-coded into a database file. This means, for this implementation, there is a maximum crossing number: 12.

## DrawPD

Perhaps one of the most exciting things I did, was learn how DrawPD works. This was a function our group had been using for a very long time. Essentially it takes a PD code and outputs a beautiful knot with smooth curves and crossings visible. It remained a mystery to most of us, but it worked nonetheless. What I found out in trying to understand this, is why the function would get stuck sometimes. It got quite unbearable at larger crossing numbers and took upwards of 2 minutes for a 14 crossing knot. I decided that it was time to look at the code and make optimizations. It turns out I was able to heuristically speed up the program almost infinitely, with a visually identical output, with the small caveat that there is a infinitessimal chance the output is not the 'perfect' 'mathematical' result. For any reasonable crossing number, however, the outputs have always been identical, but much faster. I never understood the program fully, but I have a much better sense at what is going on under the hood.

## Root files

The credit for these goes to the name in the title.
init.m is the KnotTheory` initiation package. In it provides my updated DrawPD function.