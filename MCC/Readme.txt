Minutia Cylinder-Code (MCC).

=====================================================================================================

R. Cappelli, M. Ferrara, D. Maltoni, Minutia Cylinder-Code: A New Representation and Matching Technique
for Fingerprint Recognition.  IEEE Transactions on Pattern Analysis and Machine Intelligence,  32(12) 
2128-2141 (2010).
=====================================================================================================

A local fingerprint matching algorithm. It works by constructing 3D-data structures called cylinders from the minutiae angles and position.
It supports the standard ISO/IEC 1979-2 (not the format, but the exclusive matching based on only minutiae). It is robuts to deformations and works
with smooth borders for takint into account the relative properties among the closest minutiae.


Phase I) Creates all the cylinders from the minutiae set. Each cylinder is represented by a floating or binary vector depending on the algorithm used
	in the construction. A cylinder is considered non-valid if it has not enough relation with other cylinders.

Phase II) Calculates the similarities matrix among all the valid cylinders.

Phase III) Consolidation: The algorithm implements for types of consolidation.
	- LSS: Is just the search for the best pairs of cylinders without taking account repetitions of assignment.
	- LSA: The assignment is done through the hungarian algorithm.
	- LSS-R: The similarities ar adjusted by a reinforcement algorithm based in the relation of the minutiae with other candidate minutiae. The 
		assignment is done by LSS.
	- LSA-R: Is just the combination of LSA and the reinforcement algorithm.

