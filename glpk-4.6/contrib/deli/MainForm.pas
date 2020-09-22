// ---------------------------------------------------------------------
// Copyright (C) 2003 Andrew Makhorin <mao@mai2.rcnet.ru>, Department
// for Applied Informatics, Moscow Aviation Institute, Moscow, Russia.
// All rights reserved.
//
// Author: Ivo van Baren, Operations Research specialist,
//         i.van.baren@freeler.nl
//
// This file is a part of GLPK (GNU Linear Programming Kit).
//
// GLPK is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// GLPK is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GLPK; see the file COPYING. If not, write to the Free
// Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
// --------------------------------------------------------------------



unit MainForm;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, Sample;

type
  TForm1 = class(TForm)
    Button1: TButton;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

implementation

{$R *.DFM}

procedure TForm1.Button1Click(Sender: TObject);
begin
 GLPKSolve;
end;

end.
