# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:32 2022

@author: M
"""

from math import sin, cos, sqrt, degrees, atan, pi, tan, atan2, asin
import numpy as np

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        
        Parameters
        ----------
        X, Y, Z : float : współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat : [float] : [stopnie dziesiętne] szerokość geodezyjna
        lon : [float] : [stopnie dziesiętne] długośc geodezyjna.
        h : [float] : [metry] wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")


    def deg2rad(self, kat):
        '''
        Funkcja przelicza kąt w stopniach na radiany.
        
        Parameters:
        ----------
        kat : [float/int] : kąt podany w stopniach
        
        Returns:
        -------
        katr : [float/int] : kąt przeliczony z argumentu 'kat' na radiany
        
        '''
        from math import pi
        katr = kat * pi / 180
        return(katr)


      
    def plh2XYZ(self, phi, lam, h):
        
        '''
        Funkcja przelicza współrzędne krzywoliniowe na prostokątne.
        
        Parameters:
        ----------
        
        phi : [float] : szerokość geograficzna punktu
        lam : [float] : długość geograficzna punktu 
        h : [float] : wysokość punktu      

        Returns:
        -------
        X : [float] : współrzędna prostokątna X punktu 
        Y : [float] : współrzędna prostokątna Y punktu 
        Z : [float] : współrzędna prostokątna Z punktu 
        
        '''

        a = self.a
        ecc2 = self.ecc2

        phi = self.deg2rad(phi)
        lam = self.deg2rad(lam)
        
        N =  a / (sqrt(1 - ecc2 * (sin(phi)) ** 2))
        
        X = (N + h) * cos(phi) * cos(lam)
        Y = (N + h) * cos(phi) * sin(lam)
        Z = (N * (1 - ecc2) + h) * sin(phi)
        
        return(X, Y, Z)

   

    def u2000(self, phi, lam):
        """
        Funkcja przeliczająca współrzedne prostokatne lokalne na płaszczyznie Gaussa-Krugera do ukladu PL-2000
        
        Parameters:
        -----------
        phi : [float] : szerokość geograficzna punktu
        lam : [float] : długość geograficzna punktu 
        
        Returns:
        --------
        x : [float : ]wspołrzędna x w układzie PL-2000
        y: [float] : współrzędna y w układzie PL-2000
        
        """
        
        L0=21*pi/180
        phi=self.deg2rad(phi)
        lam=self.deg2rad(lam)
        b2=(self.a**2)*(1-self.ecc2)
        ep2=((self.a**2-b2))/(b2)
        t=tan(phi)
        n2=ep2*((cos(phi))**2)
        N =  self.a / (sqrt(1 - self.ecc2 * (sin(phi)) ** 2))
        A0=1-(self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256);
        A2=(3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128));
        A4=(15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4));
        A6=(35*(self.ecc2**3))/3072;
        si=self.a*(A0*(phi)-A2*sin(2*phi)+A4*sin(4*phi)-A6*sin(6*phi));
        dl = lam-L0
        xgk=si+((dl**2)/2)*N*sin(phi)*cos(phi)*(1+((dl**2)/12)*((cos(phi))**2)*(5-t**2+9*n2+4*(n2**2))+((dl**4)/360)*((cos(phi))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
        ygk=dl*N*cos(phi)*(1+((dl**2)/6)*((cos(phi))**2)*(1-t**2+n2)+((dl**4)/120)*((cos(phi))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))
        x=xgk*0.999923
        y=ygk*0.999923+(degrees(L0)/3)*1000000+500000
        return x,y


    def u1992(self, phi, lam):
        """
        Funkcja przeliczająca współrzedne prostokatne lokalne na płaszczyznie Gaussa-Krugera do ukladu PL-2000
        
        Parameters:
        -----------
        phi : [float] : szerokość geograficzna punktu
        lam : [float] : długość geograficzna punktu 
        
        Returns:
        --------
        x : [float : ]wspołrzędna x w układzie PL-1992
        y : [float] : współrzędna y w układzie PL-1992
        
        """
        L0=19*pi/180
        phi=self.deg2rad(phi)
        lam=self.deg2rad(lam)
        b2=(self.a**2)*(1-self.ecc2)
        ep2=((self.a**2-b2))/(b2)
        t=tan(phi)
        n2=ep2*((cos(phi))**2)
        N =  self.a / (sqrt(1 - self.ecc2 * (sin(phi)) ** 2))
        A0=1-(self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256);
        A2=(3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128));
        A4=(15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4));
        A6=(35*(self.ecc2**3))/3072;
        si=self.a*(A0*(phi)-A2*sin(2*phi)+A4*sin(4*phi)-A6*sin(6*phi));
        dl = lam-L0
        xgk=si+((dl**2)/2)*N*sin(phi)*cos(phi)*(1+((dl**2)/12)*((cos(phi))**2)*(5-t**2+9*n2+4*(n2**2))+((dl**4)/360)*((cos(phi))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
        ygk=dl*N*cos(phi)*(1+((dl**2)/6)*((cos(phi))**2)*(1-t**2+n2)+((dl**4)/120)*((cos(phi))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))
        x=xgk*0.9993-5300000;
        y=ygk*0.9993+500000;
        return x,y
    
    def clean_line(self, linia, sep=','):
        """
        Funkcja służąca do odeseparowywania danych od siebie, czyszczenia lini ze spacji
        i zapisywania żądanych elementów do listy
    
        Parameters:
        ----------
        linia : linia tekstu
        sep : separator
    
        Returns:
        -------
        elementy : lista z potrzebnymi elementami
    
        """
        elementy=[]
        d=linia.split(sep)
        
        for i in d:
            if i !='':
                elementy.append(i)
        return elementy

    def wczytanie(self, dane):
        """
        Funkcja służąca do odczytania danych z pliku tekstowego i utworzenia 
        tablicy z danymi oraz słownika z danymi z nagłówka
    
        Parameters
        ----------
        dane : plik tekstowy
    
        Returns
        -------
        xyz_arr: [array of float64]: tablica z danymi
    
        """
        #otwarcie
        plik=open(dane,'r')
        xyz=[]
        #odczyt z pliku
        dane=plik.readlines()
        for linia in dane:
            if linia[0]=='3':
                a=self.clean_line(linia)
                #print(a)
                x=a[0]
                y=a[1]
                z=a[2]
                xyz.append([float(x),float(y),float(z)])
        xyz_arr=np.array(xyz)
        #print(xyz_arr)
        return xyz_arr
    
    def srednia(self, wartosci):
        """
        Funkcja liczy średnią wartość z elementów w liscie
        
        Parameters:
        ----------
        wartosci : [float] : lista wartosci
        
        Returns:
        --------
        srednia : [float] : średnia arytmetyczna elementów z listy 
        
        """
        suma = 0
        ile = 0
        for wartosc in wartosci:
            suma += wartosc
            ile += 1
        srednia = float(suma / ile)
        return(srednia)
    
    def Rneu(self, fi, l):
        """
        Funkcja, która, przyjmujac współrzedne krzywoliniowe utworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych do układu współrzędnych neu
    
        INPUT:
        ----------
        fi : [float] : wspołrzędna fi punktu początkowego układu lokalnego
        l : [float] :wspołrzędna l punktu początkowego układu lokalnego
    
        OUTPUT:
        -------
        R : [array of float64] : macierz obrotu
    
        """
        n=[(-np.sin(fi) * np.cos(l)), (-np.sin(fi) * np.sin(l)), (np.cos(fi))]
        e=[(-np.sin(l)), (np.cos(l)),  (0)]
        u=[( np.cos(fi) * np.cos(l)), ( np.cos(fi) * np.sin(l)), (np.sin(fi))]
        R=np.transpose(np.array([n,e,u]))
        return (R)
    
    def delta_neu(self, R,v):
        """
        Funckja obliczająca wektor w układzie neu
    
        Parameters:
        -----------
        R : R : [array of float64] : macierz obrotu
        v : [array of float64] : wektor w układzie XYZ
        
        Returns:
        -------
        delta_neu : [array of float64] : współrzedne topocentryczne (North (N), East (E), Up (U))
    
        """
        delta_neu=np.zeros(v.shape)
        for a in range(v.shape[0]):
            for b in range(3):
                for c in range(3):
                    delta_neu[a,c]+=v[a,b]*R[c,b]
        return (delta_neu)
    
    def odl3D(self, X1, Y1, Z1, X2, Y2, Z2) :
        """
        Funkcja do obliczenia długości na podstawie współrzędnych 3D.
        
        Parameters:
        --------------
        X1 : [float] : współrzędne X pkt 1 
        Y1 : [float] : współrzędna Y pkt 1 
        Z1 : [float] : współrzędna Z pkt 1 
        X2 : [float] : współrzędne X pkt 2 
        Y2 : [float] : współrzędna Y pkt 2 
        Z2 : [float] : współrzędna Z pkt 2 
        
        Returns:
        --------------
        odl : [float] : odległość 1-2 [m]
    
        """
        odl = sqrt( (X2 - X1)**2 + (Y2 - Y1)**2 + (Z2 - Z1)**2 )
        return(odl)
    
    def odl2D(self, X1, Y1, X2, Y2) :
        """
        Funkcja do obliczenia długości na podstawie współrzędnych 2D.
        
        Parameters:
        --------------
        X1 : [float] : współrzędne X pkt 1 
        Y1 : [float] : współrzędna Y pkt 1 
        X2 : [float] : współrzędne X pkt 2 
        Y2 : [float] : współrzędna Y pkt 2 

        
        Returns:
        --------------
        odl - [float] - odległość 1-2 [m]
    
        """
        odl = sqrt( (X2 - X1)**2 + (Y2 - Y1)**2 )
        return(odl)
    
    def azymut_elewacja(self, N, E, U):
        """   
        Funkcja wyznacza kąt azymutu i kąt elewacji 
        na podstawie współrzędnych topocentrycznych
        
        Parameters:
        -------
        N  : [float] : wpółrzedna topocentryczna N (north) [m]
        E  : [float] : wpółrzedna topocentryczna E (east) [m]
        U  : [float] : wpółrzedna topocentryczna U (up) [m] 
       
        Returns:
        -------
        Az : [float] : azymut [stopnie dziesiętne]
        el : [float] : kąt elewacji [stopnie dziesiętne]
        
        """  
        az = atan2(E, N)
        el = asin(U/sqrt(N**2 + E**2 + U**2))
        
        az=degrees(az)
        el=degrees(el)
        
        return(az, el)
    
    
    def menu_wybor(self, n):
        """
        Fukcja generująca menu do wyboru transformacji.
        """
        znak = True
        
    
        if n == "1":
            print(wspol_krzyw_lin_arr)
            np.savetxt("wsp_out.txt", wspol_krzyw_lin_arr, delimiter=',', fmt = ['%10.5f', '%10.5f', '%10.3f'], header = f"konwersja wspolrzednych geodezyjnych XYZ --> PLH\n {'fi':<10s}{'la':^10s}{'h':^10s}")
            
        elif n == "2":
            print(xyz_arr2)
            np.savetxt("wsp_out.txt", xyz_arr2, delimiter=',', fmt = ['%10.5f', '%10.5f', '%10.5f'], header = f"konwersja wspolrzednych geodezyjnych PLH --> XYZ\n {'X':<10s}{'Y':^10s}{'Z':^10s}")
        
        elif n == "3":
            print(XY2000_arr)
            np.savetxt("wsp_out.txt", XY2000_arr, delimiter=',', fmt = ['%10.5f', '%10.5f'], header = f"konwersja wspolrzednych geodezyjnych XYZ --> XY2000\n {'X':<10s}{'Y':^10s}") 
       
        elif n == "4":
            print(XY1992_arr)
            np.savetxt("wsp_out.txt", XY1992_arr, delimiter=',', fmt = ['%10.5f', '%10.5f'], header = f"konwersja wspolrzednych geodezyjnych XYZ --> XY1992\n {'X':<10s}{'Y':^10s}")  
        
        elif n == "5":
            print(neu)
            np.savetxt("wsp_out.txt", neu, delimiter=',', fmt = ['%10.5f', '%10.5f', '%10.3f'], header = f"konwersja wspolrzednych geodezyjnych XYZ --> NEU\n {'N':<10s}{'E':^10s}{'U':^10s}")
        
        elif n == "6":
            print('odległosc 3D\n', odleglosc3D)
            
        elif n == "7":
            print('odległosc 2D\n', odleglosc2D)
            
        elif n == "8":
            print('kąty azymutu\n', Azymuty)
            print('kąty elewacji\n', Elewacje)
            
        elif n == "9":
            print("""Zakonczono dzialanie programu""")
            input('')
            znak = False
        
        else:
            print("""Nie ma takiej opcji""")
            powrot()
        return znak
            
    def menu_wejscie(self):
        print ("""==============================
                MENU             
    ==============================
    1. Współrzędne fi,l,h
    2. Współrzędne X,Y Z
    3. Współrzędne X,Y 2000
    4. Współrzędne X,Y 1992
    5. NEU
    6. Odległosc 3D
    7. Odleglosc 2D
    8. Kąty azymutu i elewacji
    9. Wyjście     """)
        return
            
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # phi, lam, h = geo.xyz2plh(X, Y, Z)
    # X1, Y1, Z1 = geo.plh2XYZ(phi, lam, h)
    # x2000, y2000 = geo.u2000(phi, lam)
    # x1992, y1992 = geo.u1992(phi, lam)
    # print(phi, lam, h)
    # print(X1, Y1, Z1)
    # print(x2000, y2000)
    # print(x1992, y1992)
    xyz_arr=geo.wczytanie('wsp_inp.txt')
    #print(xyz_arr)

    
    #transformacja wspolrzednych geocentrycznych na wspolrzedne krzywoliniowe flh
    wspol_krzyw_lin=[]
    for i in xyz_arr:
        flh =geo.xyz2plh(i[0], i[1], i[2])
        wspol_krzyw_lin.append(flh)
        
    wspol_krzyw_lin_arr=np.array(wspol_krzyw_lin)
    # print('xyz --> plh\n', wspol_krzyw_lin_arr)
    
    #transformacja odwrotna
    xyz2=[]
    for i in wspol_krzyw_lin_arr:
        XYZ =geo.plh2XYZ(i[0], i[1], i[2])
        xyz2.append(XYZ)
        
    xyz_arr2=np.array(xyz2)
    # print('plh --> XYZ\n', xyz_arr2)
    
    #transformacja do ukladu PL-2000
    XY2000=[]
    for k in wspol_krzyw_lin_arr:
        xy2000= geo.u2000(k[0],k[1])
        XY2000.append(xy2000)
    
    XY2000_arr=np.array(XY2000)
    # print('układ PL-2000\n', XY2000_arr)
    
    #transformacja do ukladu PL-1992
    XY1992=[]
    for k in wspol_krzyw_lin_arr:
        xy1992= geo.u1992(k[0],k[1])
        XY1992.append(xy1992)
    
    XY1992_arr=np.array(XY1992)
    # print('układ PL-1992\n', XY1992_arr)
    
    #transformacja neu
    listaX=[]
    listaY=[]
    listaZ=[]
    for l in xyz_arr:
        listaX.append(l[0])
        X_sr=geo.srednia(listaX)
        listaY.append(l[1])
        Y_sr=geo.srednia(listaY)
        listaZ.append(l[2])
        Z_sr=geo.srednia(listaZ)
    
    i=0
    v=np.array(np.zeros((xyz_arr.shape[0],3)))
    for i in range(0,xyz_arr.shape[0]):
        v[i,0]=X_sr-xyz_arr[i,0]
        v[i,1]=Y_sr-xyz_arr[i,1]
        v[i,2]=Z_sr-xyz_arr[i,2]
        i=i+1
    [fi_sr,lam_sr,h_sr]=geo.xyz2plh(X_sr, Y_sr, Z_sr)
    R=geo.Rneu(fi_sr,lam_sr)
    neu=geo.delta_neu(R, v)
    # print('NEU\n',neu)
    
    
    #odległosc 3D
    odleglosc3D=[]
    for i in xyz_arr:
        odl=geo.odl3D(xyz_arr[0,0],xyz_arr[0,1],xyz_arr[0,2],i[0],i[1],i[2])
        odleglosc3D.append(odl)
    
    # print('odległosc 3D\n', odleglosc3D)
    
    #odległosc 2D
    odleglosc2D=[]
    for i in xyz_arr:
        odl2=geo.odl2D(xyz_arr[0,0], xyz_arr[0,1], i[0], i[1])
        odleglosc2D.append(odl2)

    # print('odległosc 2D\n', odleglosc2D)
    
    Azymuty=[]
    Elewacje=[]
    for k in neu:
        azymut, elewacja= geo.azymut_elewacja(k[0],k[1],k[2])
        Azymuty.append(azymut)
        Elewacje.append(elewacja)

    # print('kąty azymutu\n', Azymuty)
    # print('kąty elewacji\n', Elewacje)
    
    #wybór transformacji
    znak = True
    
    while znak == True:
        
        geo.menu_wejscie()
    
        n = input('Wybierz opcje: ')
        
        znak = geo.menu_wybor(n)
    
    
    