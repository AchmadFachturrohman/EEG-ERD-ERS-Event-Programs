unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, VclTee.TeeGDIPlus,
  VCLTee.TeEngine, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart,
  VCLTee.Series, Math, VCLTee.TeeSurfa, Vcl.ComCtrls, VCLTee.TeeSurfaceTool,
  VCLTee.TeeTriSurface, VCLTee.TeeTools;

type
   myArray2 = array[0..2,-100000..100000] of extended;
  TForm1 = class(TForm)
    BPFButton: TButton;
    OpenDataButton: TButton;
    Label1: TLabel;
    Edit1: TEdit;
    Edit2: TEdit;
    STFTButton: TButton;
    Chart1: TChart;
    Edit4: TEdit;
    Edit5: TEdit;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    ClearButton: TButton;
    GroupBox2: TGroupBox;
    hanning: TRadioButton;
    hamming: TRadioButton;
    Chart3: TChart;
    Chart4: TChart;
    Label8: TLabel;
    Chart2: TChart;
    TeeGDIPlus2: TTeeGDIPlus;
    TeeGDIPlus1: TTeeGDIPlus;
    OpenDialog1: TOpenDialog;
    Series1: TLineSeries;
    Series2: TLineSeries;
    Series3: TLineSeries;
    Series4: TLineSeries;
    Series5: TLineSeries;
    Series6: TLineSeries;
    Series7: TLineSeries;
    Series8: TLineSeries;
    Series9: TLineSeries;
    MAVButton: TButton;
    Chart5: TChart;
    Series10: TLineSeries;
    Series11: TLineSeries;
    Series12: TLineSeries;
    PageControl1: TPageControl;
    TabSheet1: TTabSheet;
    TabSheet2: TTabSheet;
    TabSheet3: TTabSheet;
    GroupBox1: TGroupBox;
    TabSheet5: TTabSheet;
    DFTButton: TButton;
    Chart7: TChart;
    Chart8: TChart;
    Series15: TBarSeries;
    Series16: TBarSeries;
    Series17: TBarSeries;
    Series18: TBarSeries;
    Series19: TBarSeries;
    Series20: TBarSeries;
    TabSheet6: TTabSheet;
    HPF2Button: TButton;
    Chart9: TChart;
    Chart10: TChart;
    Series21: TLineSeries;
    Series22: TLineSeries;
    SquareButton: TButton;
    Chart11: TChart;
    Series23: TLineSeries;
    TabSheet8: TTabSheet;
    ERDERSButton: TButton;
    Chart12: TChart;
    Series24: TLineSeries;
    Chart13: TChart;
    Series25: TLineSeries;
    Series13: TTriSurfaceSeries;
    ChartTool1: TRotateTool;
    Chart6: TChart;
    Series14: TSurfaceSeries;
    ScrollBar1: TScrollBar;
    Label2: TLabel;
    Label3: TLabel;
    Edit3: TEdit;
    Label4: TLabel;
    Edit6: TEdit;
    Label9: TLabel;
    Edit7: TEdit;
    Label10: TLabel;
    Edit8: TEdit;
    Label11: TLabel;
    Edit9: TEdit;
    Label12: TLabel;
    Edit10: TEdit;
    procedure OpenDataButtonClick(Sender: TObject);
    procedure proses_stft;
    procedure dft_stft;
    procedure mav;
    procedure BPFButtonClick(Sender: TObject);
    procedure ClearButtonClick(Sender: TObject);
    procedure STFTButtonClick(Sender: TObject);
    procedure MAVButtonClick(Sender: TObject);
    procedure DFTButtonClick(Sender: TObject);
    procedure HPF2ButtonClick(Sender: TObject);
    procedure SquareButtonClick(Sender: TObject);
    procedure ERDERSButtonClick(Sender: TObject);
    procedure ScrollBar1Change(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;
  jmldata, fs  :integer;
  lebar_w, geser, jml_w, jml_stft :integer;
  datadft, min, max :integer;
  mag_stft :array [0..10000,0..10000] of real;
  sinyalwindow, window  : array [-100000..100000] of real;
  freq_stft, time_stft :array [0..1000000] of extended;
  dft_re, dft_im, inputdft, sinyaldft :array [0..1000000] of extended;
  lpf2, hpf2, square, mav, mav_fw, mav_bw, erders :array [0..1000000] of extended;
  ch, lpf, hpf, autocorr, riil, imag, hasildft1, hasildft2, fw, bw :myArray2;

implementation

{$R *.dfm}
procedure TForm1.OpenDataButtonClick(Sender: TObject);
var
  i,j : integer;
  del1, del2, val : string;
  ambildata: TStringList;
  dataeeg : textfile;
begin
  Series1.Clear; Series2.Clear; Series3.Clear;
  ambildata := TStringList.Create;
  i := 0;

  if OpenDialog1.Execute = true then
  begin
    assignfile(dataeeg,OpenDialog1.FileName);

    reset(dataeeg);
    readln(dataeeg,del1);
    ambildata.Delimiter:='[';
    readln(dataeeg,del2);
    ambildata.DelimitedText := del2;
    while not EOF (dataeeg) do
    begin
      readln(dataeeg,val);
      ambildata.Delimiter := ' ';
      ambildata.DelimitedText := val;

      ch[0,i] := strtofloat(ambildata[6]);//cz
      ch[0,i] := roundto(ch[0,i],-3);
      ch[1,i] := strtofloat(ambildata[14]);//c4
      ch[1,i] := roundto(ch[1,i],-3);
      ch[2,i] := strtofloat(ambildata[18]);//c3
      ch[2,i] := roundto(ch[2,i],-3);

      jmldata := i;
      inc(i);
    end;
    closefile(dataeeg);
  end;

  fs := 500;
  Edit1.Text := inttostr(jmldata);
  Edit2.Text := inttostr(fs);

  for i := 0 to jmldata-1 do
  begin
    for j := 0 to 2 do
    begin
      Chart1.Series[j].AddXY(i/fs,ch[j,i]);
    end;
  end;
end;

procedure TForm1.BPFButtonClick(Sender: TObject);
var
  i,j :integer;
  fcL, fcH :extended;
begin
  Series4.Clear; Series5.Clear; Series6.Clear;
  Series7.Clear; Series8.Clear; Series9.Clear;

  fcL := strtofloat(Edit7.Text);
  fcH := strtofloat(Edit8.Text);
  for i := 0 to jmldata-1 do
  begin
    for j := 0 to 2 do
    begin
      lpf[j,i] := ((((8*sqr(fs))-(2*sqr(fcL)))*lpf[j,(i-1)])-(((4*sqr(fs))-((2*sqrt(2)*fcL)*fs)+(sqr(fcL)))*lpf[j,(i-2)])+(sqr(fcL)*ch[j,i])+(2*sqr(fcL)*ch[j,(i-1)])+(sqr(fcL)*ch[j,(i-2)]))/((4*sqr(fs))+((2*sqrt(2)*fcL)*fs)+(sqr(fcL)));
      hpf[j,i] := (4*lpf[j,i]-8*lpf[j,(i-1)]+4*lpf[j,(i-2)]-(2*sqr(1/fs)*sqr(fcH)-8)*hpf[j,(i-1)]-(sqr(fcH)*sqr(1/fs)-sqrt(2)*fcH*2*(1/fs)+4)*hpf[j,(i-2)])/(sqr(fcH)*sqr(1/fs)+sqrt(2)*fcH*2*(1/fs)+4);
      Chart3.Series[j].AddXY(i/fs,lpf[j,i]);
      Chart4.Series[j].AddXY(i/fs,hpf[j,i]);
    end;
  end;
end;

procedure TForm1.DFTButtonClick(Sender: TObject);
var
  i,j,k :integer;
begin
  k := 0;
  while k<=2 do
  begin
    for i := 0 to jmldata-1 do
    begin
      riil[k,i] := 0;
      imag[k,i] := 0;
      for j := 0 to jmldata-1 do
      begin
        riil[k,i] := riil[k,i] + ch[k,j]*cos(2*pi*j*i/jmldata);
        imag[k,i] := imag[k,i] - ch[k,j]*sin(2*pi*j*i/jmldata);
      end;
      hasildft1[k,i] := sqrt(sqr(riil[k,i]) + sqr(imag[k,i]))/jmldata;
    end;
    inc(k);
  end;

  k := 0;
  while k<=2 do
  begin
    for i := 0 to jmldata-1 do
    begin
      riil[k,i] := 0;
      imag[k,i] := 0;
      for j := 0 to jmldata-1 do
      begin
        riil[k,i] := riil[k,i] + hpf[k,j]*cos(2*pi*j*i/jmldata);
        imag[k,i] := imag[k,i] - hpf[k,j]*sin(2*pi*j*i/jmldata);
      end;
      hasildft2[k,i] := sqrt(sqr(riil[k,i]) + sqr(imag[k,i]))/jmldata;
    end;
    inc(k);
  end;

  for i := 0 to round(jmldata/2) do
  begin
    for j := 0 to 2 do
    begin
      Chart8.Series[j].AddXY(i*fs/jmldata,hasildft2[j,i]);
      Chart7.Series[j].AddXY(i*fs/jmldata,hasildft1[j,i]);
    end;
  end;
end;

procedure TForm1.MAVButtonClick(Sender: TObject);
var
  i, j, n: integer;
  temp :extended;
begin
  Series10.Clear; Series11.Clear; Series12.Clear;
  n:= 0;
  while n<=2 do
  begin
    for i := 0 to jmldata-1 do
    begin
      temp := 0;
      for j := 0 to 600 do
      begin
        temp := temp + sqr(hpf[n,i-j]);
      end;
      fw[n,i] := temp;
    end;
    inc(n);
  end;

  n := 0;
  while n<=2 do
  begin
    for i := 0 to jmldata-1 do
    begin
      temp := 0;
      for j := 0 to 600 do
      begin
        temp := temp + fw[n,i+j];
      end;
      bw[n,i] := temp;
    end;
    inc(n);
  end;

  for i := 0 to jmldata-1 do
  begin
    for j := 0 to 2 do
    begin
      Chart5.Series[j].AddXY(i,bw[j,i]);
    end;
  end;
end;

procedure TForm1.STFTButtonClick(Sender: TObject);
begin
  Series13.Clear; Series14.Clear;
  proses_stft;
end;

procedure TForm1.proses_stft;
var
  i, j  :integer;
begin
  lebar_w  := strtoint(Edit4.Text);
  geser    := strtoint(Edit5.Text);
  jml_stft := strtoint(Edit10.Text) - strtoint(Edit9.Text);
  jml_w    := round(jmldata/((lebar_w div 2)+geser))+1;
  Label7.Caption := 'Jumlah Window : ' + inttostr(jml_w);

  for i := 0 to jml_w-1 do
  begin
    min := i*((lebar_w div 2) + geser) - (lebar_w div 2);
    max := min + lebar_w;

    if min < 0 then
    begin
      min := 0;
    end;

    for j := 0 to jml_stft-1 do
    begin
      window[j] := 0; // inisialisasi nilai 0 untuk proses window
    end;

    for j := min to max do
    begin
      if hanning.Checked = true then  //hanning window
      begin
        window[j] := 0.5 - (0.5*cos((2*pi*(j-i*((lebar_w div 2) + geser)))/lebar_w))
      end
      else if hamming.Checked = true then   //hamming window
      begin
        window[j] := 0.54 - (0.46*cos((2*pi*(j-i*((lebar_w div 2) + geser)))/lebar_w))
      end;

      sinyalwindow[j] := bw[2,j]*window[j];
      inputdft[j-min] := sinyalwindow[j];
    end;

    datadft := max - min + 1;
    dft_stft;

    for j := 0 to round(fs/2) do
    begin
      mag_stft[i,j] := sinyaldft[j];
    end;
  end;

  for i := 0 to jml_w-1 do
  begin
    time_stft[i] := (((lebar_w div 2) + geser)*i + 22280)/fs;
  end;

  for i := 0 to jml_w-1 do
  begin
    for j := 0 to round(fs/2) do
    begin
      Series13.AddXYZ(time_stft[i],mag_stft[i,j],freq_stft[j]);
      Series14.AddXYZ(time_stft[i],mag_stft[i,j],freq_stft[j]);
    end;
  end;
end;

procedure TForm1.dft_stft;
var
  i, j :integer;
begin
  for i := 0 to jml_stft-1 do
  begin
    sinyaldft[i] := 0;
    dft_re[i] := 0;
    dft_im[i] := 0;
  end;

  for i := 0 to jml_stft-1 do
  begin
    for j := 0 to jml_stft-1 do
    begin
      dft_re[i] := dft_re[i] + inputdft[j]*cos(2*pi*j*i/jmldata);
      dft_im[i] := dft_im[i] - inputdft[j]*sin(2*pi*j*i/jmldata);
    end;
    sinyaldft[i] := sqrt(sqr(dft_re[i]) + sqr(dft_im[i]))/jmldata;
    freq_stft[i] := i*fs/datadft;
  end;
end;

procedure TForm1.HPF2ButtonClick(Sender: TObject);
var
  i :integer;
  fcL, fcH :extended;
begin
  Series21.Clear; Series22.Clear;
  fcL := strtofloat(Edit3.Text);
  fcH := strtofloat(Edit6.Text);
  for i := 0 to jmldata-1 do
  begin
    lpf2[i] := ((((8*sqr(fs))-(2*sqr(fcL)))*lpf2[i-1])-(((4*sqr(fs))-((2*sqrt(2)*fcL)*fs)+(sqr(fcL)))*lpf2[i-2])+(sqr(fcL)*ch[2,i])+(2*sqr(fcL)*ch[2,(i-1)])+(sqr(fcL)*ch[2,(i-2)]))/((4*sqr(fs))+((2*sqrt(2)*fcL)*fs)+(sqr(fcL)));
    hpf2[i] := (4*lpf2[i]-8*lpf2[i-1]+4*lpf2[i-2]-(2*sqr(1/fs)*sqr(fcH)-8)*hpf2[i-1]-(sqr(fcH)*sqr(1/fs)-sqrt(2)*fcH*2*(1/fs)+4)*hpf2[i-2])/(sqr(fcH)*sqr(1/fs)+sqrt(2)*fcH*2*(1/fs)+4);
    Series21.AddXY(i/fs,lpf2[i]);
    Series22.AddXY(i/fs,hpf2[i]);
  end;
end;

procedure TForm1.SquareButtonClick(Sender: TObject);
var
  i :integer;
begin
  Series23.Clear;
  for i := 0 to jmldata-1 do
  begin
    square[i] := sqr(hpf2[i]);
    Series23.AddXY(i/fs,square[i]);
  end;
end;

procedure TForm1.ScrollBar1Change(Sender: TObject);
begin
  mav;
end;

procedure TForm1.mav;
var
  i,j,k :integer;
  temp :extended;
begin
  series24.Clear;
  k := scrollbar1.Position;
  label2.Caption := 'Orde: ' + inttostr(k - 1);
  for i := 0 to jmldata-1 do
  begin
    temp := 0;
    for j := 0 to k-1 do
    begin
      temp := temp + square[i-j];
    end;
    mav_fw[i] := temp/k;
  end;

  for i := 0 to jmldata-1 do
  begin
    temp := 0;
    for j := 0 to k-1 do
    begin
      temp := temp + mav_fw[i+j];
    end;
    mav_bw[i] := temp/k;
  end;

  for i := 0 to jmldata-1 do
  begin
    series24.AddXY(i/fs,mav_bw[i]);
  end;
end;

procedure TForm1.ERDERSButtonClick(Sender: TObject);
var
  i :integer;
  reff :extended;
begin
  Series25.Clear;
  reff := 0;
  for i := 500 to jmldata-1 do
  begin
    reff := mav_bw[i-fs];
    if reff <> 0 then
    begin
      erders[i] := ((mav_bw[i] - reff)/reff)*100;
      Series25.AddXY(i/fs,erders[i]/10);
    end;
  end;
end;

procedure TForm1.ClearButtonClick(Sender: TObject);
begin
  Series1.Clear; Series2.Clear; Series3.Clear;
  Series4.Clear; Series5.Clear; Series6.Clear;
  Series6.Clear; Series6.Clear; Series6.Clear;
  Series7.Clear; Series8.Clear; Series9.Clear;
  Series10.Clear; Series11.Clear; Series12.Clear;
  Series13.Clear; Series14.Clear; Series15.Clear;
  Series16.Clear; Series17.Clear; Series18.Clear;
  Series19.Clear; Series20.Clear; Series21.Clear;
  Series22.Clear; Series23.Clear; Series24.Clear;
  Series25.Clear;
  Label7.Caption := ''; Label2.Caption := '';
  scrollbar1.Position := 0;
  hanning.Checked := false;
  hamming.Checked := false;
  Edit1.Clear; Edit2.Clear;
end;

end.
