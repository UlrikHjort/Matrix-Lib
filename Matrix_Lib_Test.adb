-------------------------------------------------------------------------------
--                                                                           --
--                               Matrix Lib                                  --
--                                                                           --
--                           Matrix_Lib_Test.adb                             --
--                                                                           --
--                                 MAIN                                      --
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
with Ada.Text_IO; use Ada.Text_IO;
with Matrix_Lib; use Matrix_Lib;

procedure Matrix_Lib_Test is


   M4 : constant Matrix_Lib.Matrix := (
                   (2.0,1.0,-1.0, 8.0),
                   (-3.0,-1.0,2.0,-11.0),
                   (-2.0,1.0,2.0,-3.0));



   M6 : constant Matrix_Lib.Matrix :=(
                             (1.0,3.0,3.0,1.0,0.0,0.0),
                             (1.0,4.0,3.0,0.0,1.0,0.0),
                             (1.0,3.0,4.0,0.0,0.0,1.0));


   M7 : constant Matrix_Lib.Matrix :=(
                             (1.0,3.0,3.0),
                             (1.0,4.0,3.0),
                             (1.0,3.0,4.0));




   Solution : Column(M6'Range(1));

   MI : constant Matrix := Get_Inverse(M7);
   MM : Matrix(M7'Range(1), M7'Range(2));

begin
   Put_Line("---------Test Inverse matrix:");
   Put_Line("M : ");
   Print (M7);
   Put_Line("M Inverse : ");
   Print(MI);
   Put_Line("M * Minverse : ");
   Print(M7*MI);
   Put_Line("---------End: Test Inverse matrix:");

   Put_Line("---------Test: Solve equation system:");
   Print(M4);
   Solve_Equation_System(M4, Solution);
   Put_Line("Solution: ");
   for I in Solution'Range loop
      Put_Line("X" & Integer'Image(I) & "= " & Float'Image(Solution(I)) & " ");
   end loop;
   New_Line;
   Put_Line("---------End: Test: Solve equation system:");


   Put_Line("---------Test Inverse matrix by gaussian elimination :");
   Put_Line("M : ");
   Print (M7);
   Put_Line("M Inverse : ");
   MM := Get_Inverse_By_Gaussian_Elimination(M7);
   Print(MM);
   Put_Line("M * Minverse : ");
   Print(M7*MM);
   Put_Line("---------End: Test Inverse matrix by gaussian elimination:");
end Matrix_Lib_Test;
