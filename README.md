# Octree

https://en.wikipedia.org/wiki/Octree

What is it for?
- (1) 3D spacial indexing (quick search for a point position within octree)
- (2) implicit 3D surface representation (opposite to parametric)
- (3) implicit solid geometry representation (opposite to b-rep)
- (4) geometry morphing by level sets (from this shape to this shape)
- (5) continuous function approximation on non-conformal mesh (octree is actually a non-conformal mesh) with XFEM
- (6) 3D mesh generation : create an octree with a body inside, refine octree cells near body surface and you will get non-conformal meshes inside and outside (not body-fitted yet)
- (7) conversion from a point cloud into parametric geometry

Traditional octree construction
===============================

A traditional octree construction. Every cell face (total 6) has pointers to cell neighbours.<br /><br />
Let's say you want to refine a cell (subdivide into 8 sub-cells). First check if can be refined to keep octree balanced (every two neighbour cells can have 1 level (twice size) difference. Let's say, this check is OK for new cells (after subdivision) and you make cell subdivision into 8 sub-cells.<br /><br />
The beauty of CAD programming is a false impression that everything can be done easily. You start writing the code. In this case, you need to assign proper pointers to all the old and new cells. After some testing you see that not all variants are covered and the code needs fixing. When the full code starts working as completed, again, a customer reports that in same place something is not right. You see that, again, not all variants are covered etc.<br /><br />
The conclusions is that there is a class of algorithms which seem programmable but they are actually not. They are very common in 3D.

This octree construction
=======================

No neighbour information in memory at all - all generated on the fly. Cells are fully defined by their integer coordinates of their centres. List of cells is a map of three integer coordinates into an instance of cell class. If you wish to refine a cell, you generate integer coordinates of its neighbours and add them to the map. That's all.<br /><br /> 
This makes the code very simple, reliable and saves memory. Downside : yes, it must be slower but not very much.

More detail
===========
Actually octree has two maps, one for cells and one for nodes. Nodes are 8 nodes of each cell around its centre.<br /><br />
Nodes are needed to build a good continuous approximation of a function (e.g. level-set function) which can be done with conforming finite elements (https://github.com/AndreyKoudr/FiniteElements).<br /><br />
When you add a cell to an octree, 8 new nodes are inserted into the node map, increasing their reference counts. When a cell is excluded, its nodal ref count is  decreased and a node is physically removed from total node map if only its ref count reaches zero.

Octree implementation
=====================

  <B>Background</B><br />

  Background is a set of largest cells obtained by uniform division of the whole cuboid region.
  Background cells are those of level 0. Every cell of level 0 can be subdivided (refined)
to make a hierarchy of cells inside a background cell. 

  <B>Search</B><br />

  A search of every 3D point inside octree is two-step :
  (1) which background cell? 
  (2) starting from this background cell, search within a hierarchy
Both operations are fast. This is the way an octree provides a kind of indexing for 3D points.<br />
  Use <I>findCellCentreAndCheck()</I> to find a cell for a point position.

  <B>Cells and nodes</B><br />

  Octree is a collection of cells (<I>OCELLS cells_;</I>) and nodes (ONODES nodes_;) both wholly 
defined by their single integer coordinates (cell centre for cells) <I>IPosition</I>. 
Position of cell also uniquely defines its level (number of cell subdivisions from background 
cell) (use <I>level()</I> for that) and its size (<I>intCellSize()</I>).<br />
  
  Each cell and every node carries data in OCELL_DATA and ONODE_DATA type variables. 

  <B>Cell neighbours</B><br />

  No information about neighbours is in memory inside a cell : all info is generated on the 
fly by search of neighbour cells in <I>OCELLS</I> map. This makes the code very simple, reliable and 
saves lots of memory.

  <B>Refine/derefine cells</B><br />

  Cells can be refined and de-refined. Max refinement level (for the the smallest posiible cell) is 
defined by <I>maxLevel</I> in the background. If you specify maxLevel as 20, it means that the smallest 
cell may have size of 

  background_cell_size / 2^20

that is, ~0.001mm for background size 1m.

  Call <I>refineCell()</I> to subidvide a cell into 8 sub-cells; it is not ALWAYS possible because all newly generated
cells and their neighbours must satisfy the <B>balanced</B> ocree condition : difference in levels of any two 
neigbour cells must not be higher than 1.<br />
  Call <I>derefineCell()</I> on a parent cell to delete all its sub-cells. <I>isLeaf()</I> function tells if
a cell has children or not (is a leaf). DE-refinement may not be successful as well due to a failing 
balanced condition.

Simple improvements
===================

Maps in 

	#define OCELLS std::map<IPosition,OCell<OCELL_DATA>,IPosCompare>
	#define ONODES std::map<IPosition,ONode<ONODE_DATA>,IPosCompare>
	
can be replaced by Google btree::maps (not Google maps) from cpp-btree-1.0.1. It will save ~40% of memory. I once created an octree with 250 million cells on a 16G laptop in 6 mins.

Files and classes
=================

- <B>Types.h</B>. Some useful types and macros
- <B>Vector.h</B>. These are two basic classes for 3D vectors with operators which make
possible to use simple operations on vectors like v3 = v1 + v2 etc. (https://github.com/AndreyKoudr/3DVector). Only first class without SIMD with three integet components is used here :<br />

		/** Very important class to define integer positions of cells and nodes within octree. */
		#define IPosition TVector<LINT>
	
and one less important with floating point components :

	/** Octree real type */
	#define OREAL double

	/** Octree vector type */
	#define OVECTOR TVector<OREAL>
	
- <B>Strings.h</B>. Basic string operations.

- <B>Triangles.h</B>. Auxiliary class (https://github.com/AndreyKoudr/3DTriangles) to convert octree faces into triangles with <I>void cellsToTriangles()</I>, save in an STL file for display. Nothing more.

- <B>ONode.h</B>. ONode class. A node for a finite-element linear conforming approximation within octree. Can also carry data of type T.
When you add a cell to an octree, 8 new nodes are inserted into the node map, increasing their reference counts.
When a cell is excluded, its nodal ref counts are decreased and nodes are physically removed from total node map is only ref count reaches zero.

- <B>OCell.h</B>. OCell class. Cell (octant) identified by its integer position in a map. It can also carry a data of type T.

- <B>OBackground.h</B>. Background parameters for a balanced octree. Background is a set of largest cells obtained by uniform division of the whole cuboid region.
Background cells are cells of level 0. Every cell of level 0 can be subdivided (refined) to make a hierarchy of cells inside a background cell.<br />

  A search of every 3D point inside octree is two-step :<br />
  (1) which background cell?<br />
  (2) starting from this background cell, search within a hierarchy<br />
Both operations are fast. This is the way an octree provides a kind of indexing for 3D points.

- <B>OOctree.h</B>. Octree implementation with abundant comments.

Projects
========
- project <I>STLToPointCloud</I> is an auxiliary project to generate sample point clouds from STL files. 
It has nothing in common with octrees.
- project <I>PointCloudToSTL</I> is used to demonstrate how octree is working. Its output is not always perfect;
it would require more time for better results. You must remember that the <B>Octree</B> is the main issue.<br />
See comments inside the code.<br />
<B>Please unzip the file shuttle.0.0005 in PointClouds directory!</B> if you wish to run Debug with command-line parameters in Visual Studio. This file actually has resolution 0.001, not 0.0005. I had to replace it due to Git file upload limitation.

What is good and what is bad
============================
<B>OOctree</B> class is reliable and can be used for many different purposes.<br />
Do not expect perfect watertight output in <I>PointCloudToSTL</I> project. Read explanations in <I>PointCloudToSTL.cpp</I>.




