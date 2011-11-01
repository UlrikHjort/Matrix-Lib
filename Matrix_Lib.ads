-------------------------------------------------------------------------------
--                                                                           --
--                               Matrix Lib                                  --
--                                                                           --
--                             Matrix_Lib.ads                                --
--                                                                           --
--                                 SPEC                                      --
--                                                                           --
--                   Copyright (C) 1996 Ulrik Hørlyk Hjort                   --
--                                                                           --
--  Matrix Lib is free software;  you can  redistribute it                   --
--  and/or modify it under terms of the  GNU General Public License          --
--  as published  by the Free Software  Foundation;  either version 2,       --
--  or (at your option) any later version.                                   --
--  Matrix Lib is distributed in the hope that it will be                    --
--  useful, but WITHOUT ANY WARRANTY;  without even the  implied warranty    --
--  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  --
--  See the GNU General Public License for  more details.                    --
--  You should have  received  a copy of the GNU General                     --
--  Public License  distributed with Yolk.  If not, write  to  the  Free     --
--  Software Foundation,  51  Franklin  Street,  Fifth  Floor, Boston,       --
--  MA 02110 - 1301, USA.                                                    --
--                                                                           --
-------------------------------------------------------------------------------
package Matrix_Lib is
   type Matrix is array (Natural range <>, Natural range <>) of Float;
   type Column is array (Natural range <>) of Float;

   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   function "*" (Left, Right : Matrix) return Matrix;

   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   function "+" (Left, Right : Matrix) return Matrix;

   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   function "-" (Left, Right : Matrix) return Matrix;

   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   procedure Print(M : Matrix);


   ----------------------------------------------------------------------------
   --
   --  Solve the equation system AX=0
   --
   --  X is returned in Column
   ----------------------------------------------------------------------------
   procedure Solve_Equation_System(A : in Matrix; Solution : out Column);


   ----------------------------------------------------------------------------
   --
   -- Calculate the determinant for the n*n matrix A by the Leibniz formula
   -- as the sum of all permutations P over the set S in A'Range(1):
   --
   -- Det = SUM(parity(P) PROD(AiPi))
   --        P in S       i = 1 .. n
   --
   ----------------------------------------------------------------------------
   function Determinant(A : in Matrix) return Float;

  ----------------------------------------------------------------------------
  --
  -- Gaussian eliminination of the matrix A
  --
  ----------------------------------------------------------------------------
   procedure Gaussian_Elimination(A : in out Matrix);


  -------------------------------------------------------------------------------
  -- Inverse matrix I of A calculated by Gaussian elimination by solving the
  -- equation system:
  --
  --     AX = Bi, where Bi is the i row in the identity matrix of same size as A
  --              and B is size (A_row, 1)
  --
  -- For each solution Bi, Ii = Bi
  --
  --
  -- Returns the inverse matrix to A
  -------------------------------------------------------------------------------
   function Get_Inverse_By_Gaussian_Elimination(A : Matrix) return Matrix;

  -------------------------------------------------------------------------------
  --
  -- Calculate the inverse of matrix A
  --
  -------------------------------------------------------------------------------
   function Get_Inverse(A : Matrix) return Matrix;

end Matrix_Lib;
