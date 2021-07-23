# Octree

Octrees
=======

https://en.wikipedia.org/wiki/Octree

What is it for?
- (1) 3D spacial indexing (quick search for a point position within octree)
- (2) implicit 3D surface representation (opposite to parametric)
- (3) implicit solid geometry representation (opposite to b-rep)
- (4) geometry morphing by level sets (from this shape to this shape)
- (5) continuous function approximation on non-conformal mesh (octree is actually a non-conformal mesh) with XFEM
- (6) 3D mesh generation : create an octree with a body inside, refine octree cells near body surface and you will get non-conformal meshes inside and outside (not body-fitted yet)

Traditional octree construction
===============================

A traditional octree construction. Every cell face (total 6) has pointers to cell neighbours.<br />
Let's say you want to refine a cell (subdivide into 8 sub-cells). First check if can be refined to keep octree balanced (every two neighbour cells can have 1 level (twice size) difference. Let's say, this check is OK for new cells (after subdivision) and you make cell subdivision into 8 sub-cells.<br />
The beauty of CAD programming is an false impression that everything can be done easily. You start writing the code. In this case, you need to assign proper pointers to all the old and new cells. After some testing you see that not all variants are covered and the code needs fixing. When the full code starts working as completed, again, a customer reports that in same place something is not right. You see that, again, not all variants are covered etc.<br />
The conclusions is that there is a class of algorithms which seem programmable but they are actually not. They are very common in 3D.

This octree contruction
=======================

No neighbour information in memory at all - all generated on the fly. Cells are fully defined by their integer coordinates of their centres. List of cells is a map of three integer coordinates into cell class. If you wish to refine a cell, you generate integer coordinates of its neighbours and add them to the map.
This makes the code very simple, reliable and saves memory. Downside : yes, it must be slower but not very much.

More detail
===========
Actually octree has two maps, one for cells and one for nodes. Nodes are 8 nodes of each cell around its centre. 
Nodes are needed to build a good continuous approximation of a function (e.g. level-set function) which can be done with conforming finite elements (https://github.com/AndreyKoudr/FiniteElements).
When you add a cell to an octree, 8 new nodes are inserted into the node map, increasing their reference counts. When a cell is excluded, its nodal ref counts are decreased and a node is physically removed from total node map is only its ref count reaches zero.

