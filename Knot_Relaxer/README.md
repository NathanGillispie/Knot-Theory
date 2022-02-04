## Credits

Written by Jackson Powers and Nathan Gillispie in 2020 for their CPS project.

## Instructions
Evaluate the expression under the "Evaluate Me!" section.
In the main menu, you should see three options to select a knot. 
Click import points to select one of the included mtx files.
3_1a is identical to 3_1b but opening the files shows different ways of showing the 
data.

In the draw knot interface, only two beads exist per z-plane. Every time you click,
a point is drawn to the x and y coordinates of the click with the same z coordinate 
as the previous bead. Then another bead is drawn up/down to the incremented z value.
Left click increments the z value up, right click increments the z value then mult-
iplies it by -1. This means each time you click you are either going entirely above
or below the previous beads. This allows for any given alternating knot to be drawn.
Clicking save will return the drawn coordinates to the main menu. This module can
also be used independently.

Get knot data uses Mathematica's built in knot data. the points value is the number
of points used to interpolte the function used by mathematica. It is the user's 
responsibility to make sure the number of points is sufficient to describe the knot.

## Notes

There are no bead addition/deletion procedures implemented and the algorithm used 
for collision detection can be improved greatly. Thus, the simplification is slow
and knots will get stuck sometimes. 