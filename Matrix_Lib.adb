-------------------------------------------------------------------------------
--                                                                           --
--                               Matrix Lib                                  --
--                                                                           --
--                             Matrix_Lib.adb                                --
--                                                                           --
--                                 BODY                                      --
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
with Ada.Float_Text_IO; use Ada.Float_Text_IO;
with Number_Theory_Tools; use Number_Theory_Tools;
with Permutations_Generic;

package body Matrix_Lib is

    package Permutations_Package is new Permutations_Generic(Natural); use Permutations_Package;
   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   function "+" (Left, Right : Matrix) return Matrix is
     RetVal : Matrix(Left'Range(1), Left'Range(2)) := (others =>(others => 0.0));

   begin
      for I in Left'Range(1) loop
         for J in Left'Range(2) loop
            RetVal(I,J) :=Left(I,J) + Right(I,J);
          end loop;
      end loop;
      return RetVal;
   end "+";

   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   function "-" (Left, Right : Matrix) return Matrix is
     RetVal : Matrix(Left'Range(1), Left'Range(2)) := (others =>(others => 0.0));

   begin
      for I in Left'Range(1) loop
         for J in Left'Range(2) loop
            RetVal(I,J) :=Left(I,J) - Right(I,J);
          end loop;
      end loop;
      return RetVal;
   end "-";

   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   function "*" (Left, Right : Matrix) return Matrix is
      RetVal : Matrix(Left'Range(1), Right'Range(2)) := (others =>(others => 0.0));
   begin
      if Left'Length(2) /= Right'Length(1) then
         raise Constraint_Error;
      end if;

      for I in Left'Range(1) loop
         for J in Right'Range(2) loop
            for K in Left'Range(2) loop
               RetVal(I,J) := RetVal(I,J) + Left(I, K)*Right(K, J);
            end loop;
         end loop;
      end loop;
      return RetVal;
   end "*";


   ----------------------------------------------------------------------------
   --
   ----------------------------------------------------------------------------
   procedure Print(M : Matrix) is
      begin
      for I in M'Range(1) loop
         for J in M'Range(2) loop
            Put(Item => M(I,J), Aft => 1 ,EXP => 0); Put(" ");
         end loop;
         New_Line;
      end loop;
   New_Line;
  end Print;


  ----------------------------------------------------------------------------
  --
  ----------------------------------------------------------------------------
  procedure Forward_Substitution(A : in out Matrix)  is
     M    : constant Natural := A'Last(1);
     N    : constant Natural := A'Last(2);
     Maxi : Natural;
     Element : Float;

  begin
     for I in A'Range(1) loop
        Maxi := I;

        for J in (I + 1) .. M loop
           if A(J,I) > A(Maxi,I) then
                 Maxi := J;
              end if;
        end loop;

        -- Interchange rows:
        for J in A'Range(2) loop
           Element   := A(Maxi,J);
           A(Maxi,J) := A(I,J);
           A(I,J)    := Element;
        end loop;

        for J in reverse I .. N loop
           for K in (I + 1) .. M loop
             A(K,J) := A(K,J) - (A(K,I)/A(I,I) * A(I,J));
           end loop;
        end loop;
     end loop;
  end Forward_Substitution;

  ----------------------------------------------------------------------------
  --
  ----------------------------------------------------------------------------
  procedure Reverse_Elimination(A : in out Matrix)  is
     N : constant Natural := A'Last(2);

  begin
       for I in reverse A'Range(1) loop

        A(I,N) := A(I,N) / A(I,I);
        A(I,I) := 1.0;
        for J in reverse 0 .. (I-1) loop
           A(J,N) := A(J,N) - (A(J,I) * A(I,N));
           A(J,I) := 0.0;
        end loop;
     end loop;
  end Reverse_Elimination;


  ----------------------------------------------------------------------------
  --
  -- Gaussian eliminination of the matrix A
  --
  ----------------------------------------------------------------------------
  procedure Gaussian_Elimination(A : in out Matrix)  is

  begin
      Forward_Substitution(A);
      Reverse_Elimination(A);
  end Gaussian_Elimination;


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
  function Get_Inverse_By_Gaussian_Elimination(A : Matrix) return Matrix is
     Inverse : Matrix(A'Range(1), A'Range(2)) := (others => (others => 0.0));
     M : Matrix(A'Range(1), A'First(2) .. A'Last(2) + 1) := (others => (others => 0.0));

  begin
     for K in A'Range(2) loop

        -- Copy A into M:
        for I in A'Range(1) loop
           for J in A'Range(2) loop
              M(I,J) := A(I,J);
           end loop;
           M(I,M'Last(2)) := 0.0;
        end loop;

        -- Create next unit column:
        M(K,M'Last(2)) := 1.0;
        Gaussian_Elimination(M);

        -- Insert next column in inverse matrix:
        for L in A'Range(1) loop
           Inverse(L,K) := M(L,M'Last(2));
        end loop;
     end loop;

     return Inverse;
  end Get_Inverse_By_Gaussian_Elimination;


  ----------------------------------------------------------------------------
  --
  --  Solve the equation system AX=0
  --
  --  X is returned in Column
  ----------------------------------------------------------------------------
  procedure Solve_Equation_System(A : in Matrix; Solution : out Column) is
     A_Copy : Matrix := A;

  begin
     Gaussian_Elimination(A_Copy);
     for I in Solution'Range loop
        Solution(I) := A_Copy(I,A_Copy'Last(2));
     end loop;
  end Solve_Equation_System;




  ----------------------------------------------------------------------------
  --
  -- Calculate the determinant for the n*n matrix A by the Leibniz formula
  -- as the sum of all permutations P over the set S in A'Range(1):
  --
  -- Det = SUM(parity(P) PROD(AiPi))
  --        P in S       i = 1 .. n
  --
  ----------------------------------------------------------------------------
  function Determinant(A : in Matrix) return Float is
     Permutations_List      : Permutations_T(0 .. Natural(Factorial(A'Length(1))), A'Range(1));
     Number_Of_Permutations : Natural := 0;
     Initial                : List_T(A'Range(1));
     Permutation            : List_T(A'Range(1));
     Sign                   : Float;
     Product                : Float := 1.0;
     Sum                    : Float := 0.0;

  begin
     for I in A'Range(1) loop
        Initial(I) := I;
        Permutation(I) := I;
     end loop;

     Permute(Permutation,0, Permutations_List,Number_Of_Permutations);

     for I in Permutations_List'First .. Permutations_List'Last-1 loop
        Permutation := Get_Permutation(I, Permutations_List);
        Sign := Float(Parity(Initial,Permutation));
        Product := 1.0 * sign;
        for J in Permutation'Range loop
          Product := Product * A(J,Permutation(J));
        end loop;
        Sum := Sum + Product;
     end loop;

     return Sum;
  end Determinant;

  -------------------------------------------------------------------------------
  --
  -- Returns the cofactor Crow,column
  --
  -------------------------------------------------------------------------------
  function Get_Minor(A : Matrix; Row : Natural; Column : Natural) return Matrix is
     Minor         : Matrix(A'First(1) .. A'Last(1)-1, A'First(2) .. A'Last(2)-1) := (others =>(others =>0.0));
     Row_Number    : Natural := 0;
     Column_Number : Natural := 0;

  begin
     for I in A'Range(1) loop
        if I /= Row then
           Column_Number := 0;
           for J in A'Range(1) loop
              if J /= Column then
                 Minor(Row_Number,Column_Number) := A(I,J);
                 Column_Number := Column_Number + 1;
              end if;
           end loop;
             Row_Number := Row_Number + 1;
        end if;
     end loop;
     return Minor;
   end Get_Minor;


  -------------------------------------------------------------------------------
  --
  -- Calculate the inverse of matrix A
  --
  -------------------------------------------------------------------------------
   function Get_Inverse(A : Matrix) return Matrix is
      Inverse : Matrix(A'Range(1), A'Range(2)) := (others =>(others =>0.0));
     Minor : Matrix(A'First(1) .. A'Last(1)-1, A'First(2) .. A'Last(2)-1) := (others =>(others =>0.0));
     Det : Float := Determinant(A);

   begin
      if Det /= 0.0 then
         Det := 1.0 / Det;
        for J in A'Range(1) loop
           for I in A'Range(1) loop
              Minor := Get_Minor(A,J,I);
              Inverse(I,J) := Det * Determinant(Minor);
              if (I+J) mod 2 = 1 then
                 Inverse(I,J) := - Inverse(I,J);
              end if;
           end loop;
        end loop;
      end if;
      return Inverse;
  end Get_Inverse;
end Matrix_Lib;
