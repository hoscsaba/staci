#ifndef CSOMOPONT_H
#define CSOMOPONT_H

class Csomopont
{
public:
    /// Bemen� �gak nyilv�ntart�sa
    vector<int> ag_be;
    /// Kimen� �gak nyilv�ntart�sa
    vector<int> ag_ki;

private:
    /// A csomopont neve
    string nev;
    /// HEAD = p[Pa]/ro/g
    double p_head;
    /// Suruseg
    double ro;
    /// Geodetikus magassag
    double h;
    /// Fogyasztas
    double fogy;
    /// Betaplalt koncentracio
    double cl_be;
    /// Csomoponti atlagkoncentracio
    double konc_atlag;
    /// Csomoponti atloagos vizkor
    double tt; // travel time
    

public:
    /// Konstruktor
    Csomopont(const string nev,
              const double h,
              const double fogy,
              const double cl_be,
              const double pressure,
              const double ro,
              const double tt);
    /// Copy Konstruktor
    Csomopont(const Csomopont &csp);
    /// Destruktor
    ~Csomopont();
    // Info
    string Info();
    /// Nyomas beallitasa
    void Set_p(double x)
    {
        p_head = x;
    }
    /// Fogyasztas erteke
    double Get_fogy()
    {
        return fogy;
    }
    /// Fogyaszt�s �rt�ke
    void Set_fogy(double a_fogy)
    {
        fogy = a_fogy;
    }
    /// Az elem neve
    string Get_nev()
    {
        return nev;
    }
    /// Get head (p_head)
    double Get_p()
    {
        return p_head;
    }
    /// Geodetikus magassag erteke
    double Get_h()
    {
        return h;
    }

    /// Inicializ�ci�
    void Ini(int mode, double value);
    /// Double �rt�kek be�ll�t�sa
    void Set_dprop(string mit, double value);
    /// Double �rt�kek lek�r�se
    double Get_dprop(string mit);
    

};
#endif