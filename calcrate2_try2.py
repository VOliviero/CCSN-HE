#!/usr/bin/env python
"""
    Calcrate is a utility to access instrument response functions and fluxes.

    Calcrate can produce the following output-types (option -m) :
    -m energy : print table vs E of Aeff, Flux, efficiencies, and rates
    -m rates  : summary of rates per flavor and channel.
    -m mrf    : print a table of the total signal, background and model
                rejection factor (mrf) for each declination, and cut.
    -m opt    : ottimizza il cono e il cut logE per il miglior MRF

examples:

    calcrate2.py -m rates -F 1e-4*(E**-2) -d -90,-80,...,60

        prints rates for a pointsource with a flux of 1e-4/E^2 GeV-1 s-1 m-2 at 
        declinations -90 to 60 with steps of 10

    calcrate2.py -m mrf -d -45.6 -g -F "1.46e-10*(E/1e3)**(-1.32) * exp(-E/1.4e4)" --extension 0.8 --shape 0 -a 0.5,1.0,...,5.0

        input gamma ray spectrum (vela-x) model, and use a disk-extension of 0.8 degrees. The search cone runs
        from 0.5 to 5 degrees.

options:

 -f, irf             : RDF/ingredients ROOT files                         default=$AADIR/data/aanet_detresponse.cut2.v1.0.0.root
 -m, output_type     : output-type (see above)                            default=energy,rates
 -F, fluxexpr        : neutrino flux expression [1/(GeV m2 s sr)]         default=(1e-4*(E**-2))
 -g, gammaflux       : flux expression refers to gamma-rays               default=False
 -T, livetimes       : livetime (years)                                   default=1,
 -t, mcfiletypes     : Neutrino flavour-interactions to load IRF for      default=numuCC,anumuCC,nueCC,anueCC,muonMuon
 -n, flavors         : neutrino flavours                                  default=numu,anumu,nue,anue
 -c, channels        : observational channels                             default=track,shower
 -p, point           : point source mode (flux in /GeVm2s)                default=True
 -d, declinations    : declination(range) for point-mode (degrees)        default=0,
 -a, cones           : search cone radius(range) for point-mode (degrees) default=1.0,
 -e, logE_thresholds : minimal log of reconstructed energy (range)        default=1.0,
 -w, extension       : source width (degrees)                             default=0.0
 -u, shape           : source shape (if width>0); 0=disk, 1=gaus          default=0   
 -z, conflevel       : confidence level for limits                        default=0.90
 -s, significance    : significance (p-value) level for discovery         default=0.001
 -y, output_style    : output style (str,html,latex)                      default=str
 -P, fluxfile        : flux data file path                                [default: ./file_flux.txt]
 -W, flag_flux       : flux mode: Expr or Graph                           default=Expr
 -h                  : help message and exit                              default=False
 -i                  : Python interative mode (prompt when done)          default=False
"""

import ROOT, aa, array
import sys
from ana import search
from ROOT.defs import track, shower
from ROOT.defs import muon, numu, anumu, nue, anue, anumuCC, numuCC, nueCC, anueCC, muonMuon
from ROOT import flavor_name, channel_name
from ana.search.scripts.limits import *
import itertools
from math import log10

sec_per_year = 31557600

def bail(*argv):
    print("ERROR", *argv)
    sys.exit()

def my_graph_flux(flux_file_path):
    """Load flux data from a file and return a ROOT TGraph"""
    try:
        f = open(flux_file_path, "r")
    except Exception as e:
        bail("Impossibile aprire il file di flux:", flux_file_path, e)
    logE, flux = array.array('d'), array.array('d')
    for line in f:
        if line.startswith("#"):
            continue
        energy, flux_value = map(float, line.split())
        logE.append(log10(energy))  # log(energy)
        flux.append(flux_value)
    f.close()
    graph_flux = ROOT.TGraph(len(logE), logE, flux)
    return graph_flux

    
class CalcRate:
    
    def __init__(self, **kwargs):
        """ Create the Calcrate object """
        
        self.irf             = aa.aadir + "/ana/search/irfs/20223161526_DETECTORESPONSE_zen_CutLevel2_stoombdt_b40.root"
        self.mcfiletypes     = [numuCC, anumuCC, nueCC, anueCC, muonMuon]
        self.channels        = [track, shower]
        self.flavors         = [numu, anumu, nue, anue, muon]
        self.fluxexpr        = "(1e-4*(E**-2))"
        self.gammaflux       = False 
        self.livetimes       = [1.0]
        self.point           = True
        self.declinations    = [-60]
        self.cones           = [1.0]
        self.logE_thresholds = [1.0]
        self.extension       = 0
        self.shape           = ""
        self.conflevel       = 0.90
        self.significance    = 5.733e-7  # 5 sigma
        self.output_type     = ['rates']
        self.output_style    = 'str'
        self.output          = lambda table: print(table)
        
        # Impostiamo il tipo di flux in base alle opzioni passate:
        # - W (flux_mode) e -P (flux_file)
        
        
        self.flag_flux = kwargs.get("flag_flux", "Expr")  # Controllo che il valore sia letto
        self.fluxfile = kwargs.get("fluxfile", "./file_flux.txt")  # Deve corrispondere all'opzione -P!
        if "opt" in kwargs:
            self.output_type.append("opt")

        print(f"Using {self.flag_flux} for flux")  # Debug
        print(f"Received flux file path: {self.fluxfile}")  # Debug

        
        if self.flag_flux == "Expr":
            print("Using power-law expression for flux:", self.fluxexpr)
            # Qui si utilizzerebbe self.fluxexpr per il calcolo
        elif self.flag_flux == "Graph":
            print("Using graph flux")
            self.fluxfile = kwargs.get("fluxfile", "./file_flux.txt")
            print("Flux file path:", self.fluxfile)
            self.fluxgraph = my_graph_flux(self.fluxfile)
   




        for k,v in kwargs.items() :
            if k in ('channels','flavors','mcfiletypes') :
                if type(v) == list:
                    v = [ getattr( ROOT.defs,x) if type(x) == str else x for x in v ]
                else :
                    v = [  getattr( ROOT.defs,v) ]


                print (k,v,type)

            if hasattr( self, k ):
                #if type( getattr(self,k) ) == list and type(v) != list : v = [v] # list stays a list
                setattr( self, k, v )
            else : 
                if k in ['h','i'] : pass # fine
                else :
                    bail("unknown setting for calcrate", k, v )

        print(self)

        self.init()


    def __str__ ( self ) :
        T = ROOT.Table( "property","value")
        for k,v in self.__dict__.items () : T.add( k,v )
        return str(T)


    def init(self, detres = None ):
        
        if detres:
            self.detres = detres
        else:
            self.detres = ROOT.FullDetResponse( "detres", 
                                                self.irf, 
                                                self.mcfiletypes,
                                                self.channels,
                                                self.extension, self.shape=="gaus" )

        if self.gammaflux :
            self.flux = ROOT.NuFromGammaFlux( self.fluxexpr.replace("E","x[0]") )
        elif self.flag_flux == "Graph":
        # Se il flusso è un grafico (TGraph), usiamo quello caricato dal file
            self.flux = ROOT.GraphFlux( self.fluxgraph)
        else :
            self.flux = ROOT.NuExprFlux( 1, self.fluxexpr.replace("E","x") )
            

    def ranges(self) :
        return list( itertools.product( self.livetimes, self.declinations, self.cones, self.logE_thresholds ) )


    def call( self, func ) :
        for x in self.ranges() : 
            func( *x )

    def doit( self ) :
        for t in self.output_type :
            if t == "energy" : self.call( self.energy_table )
            if t == "rates"  : self.call( self.rates_table  )
            if t == "opt": self.print_optimized_mrf_table()  # Chiamata alla funzione per l'ottimizzazione
            if t == "mrf"    : 
                for chan in self.channels : self.mrf_table( chan )


    def get_sb_flav( self, flav, chan, livetime, declination, cone, logE_threshold ) :
        "signal and background rate per flavor"
        D = self.detres.get(flav,chan)
        s,b,_ = D.compute_rates(self.flux, livetime * sec_per_year, self.point, declination*pi/180 , cone*pi/180, logE_threshold )
        #print(f"get_sb_flav: s={s}, b={b} for cone={cone}, declination={declination}, logE_threshold={logE_threshold}")

        return s,b
    
    '''def get_sb_mod( self, flav, chan, livetime, declination, cone, logE_threshold ) :
        "signal and background rate per flavor"
        D = self.detres.get(flav,chan)
        s,b,_ = D.compute_rates(self.flux, livetime * sec_per_year, self.point, declination*pi/180 , cone*pi/180, logE_threshold )
        print(f"get_sb_mod: s={s}, b={b} for cone={cone}, declination={declination}, logE_threshold={logE_threshold}")

        lim = mean_limit( b , self.conflevel,False )
        mrf = lim/s
        mu_lds=get_mu_lds(b,self.significance,0.5)
        mdp=mu_lds/s
        return s,b,lim,mrf,mu_lds,mdp'''


    def get_sb(self, chan, livetime, declination, cone, logE_threshold ) :
        "get total signal and background rate, integrated over flavors"
        L =  [ self.get_sb_flav( flav, chan, livetime, declination, cone, logE_threshold ) for flav in self.flavors ]
        return list( map( sum, zip(*L)))


    def rates_table( self, livetime, declination, cone, logE_threshold ) :

        T = ROOT.Table("channel","flavor","signal","atm-background")
        T.title = f" summary for Tlive={livetime} yr, decl={declination}\n cone={cone} deg, LogEmin={logE_threshold}"

        for chan in self.channels :
            SS,BB = 0,0
            for flav in self.flavors :  
                S,B = self.get_sb_flav( flav, chan, livetime, declination, cone, logE_threshold )
                SS,BB = SS+S, BB+B
                T.add(  channel_name( chan ), flavor_name(flav), S, B )
            T.add(channel_name( chan ), "total" , SS, BB )

        self.output(T)
        return T


    def energy_table( self, livetime, declination, cone, logE_threshold ) :

        "Print an extensive table for every channel and flavor."

        print( self.channels )
        print( self.flavors )

        for chan in self.channels :
            for flav in self.flavors :
                T = self.detres.get(flav,chan).compute_rates_table( self.flux, livetime * sec_per_year, 
                                                                    self.point, declination * pi/180, cone * pi/180, logE_threshold )
                self.output(T)

        return T


    def mrf_table( self, chan ) :

        T = ROOT.Table("year","decl","cone","minloge","total-sig","total-bg","exp.limit","MRF","n limit discovery","MDP")
        T.title = "  summary for" + ROOT.channel_name( chan ) + " channel."

        for liv,dec,con,thr in self.ranges() :

            sig,bg = self.get_sb( chan, liv,dec,con,thr )
    
            lim = mean_limit( bg, self.conflevel, False )
            dis = get_mu_lds( bg, self.significance, 0.5 )
            T.add( liv, dec, con, thr, sig, bg, lim, lim/sig ,dis, dis/sig )

        self.output(T)
        return T

    def print_optimized_mrf_table(self):
    
        for chan in self.channels:  # Itera sui canali definiti
            for flav in self.flavors :
                table = ROOT.Table("livetime", "decl", "cone(opt)", "minlogE(opt)", "total-sig", "total-bg", "mean limit", "MRF", "n limit discovery", "MDP")
                
                # Se outfile non è 'interact', scrivi su file
                if self.output_style != "interact":
                    out_txt = open(self.output_style + "_optimize_MRF.dat", "w")

                def score(cone, logE_threshold):
                    # Usa get_sb_flav per ottenere s e b per il canale/flavor specifico
                    s, b = self.get_sb_flav(flav, chan, self.livetimes[0], self.declinations[0], cone, logE_threshold)
                    '''if s == 0:
                        print(f"Warning: s=0 per cone={cone}, minloge={minloge}")
                        return 1e12'''
                    #print(f"Direct call: flav={flav}, chan={chan}, livetime={livetime}, declination={declination}, cone={cone}, minloge={logE_threshold}")
                    #print(f"Inside minimize: flav={flav}, chan={chan}, livetime={self.livetimes[0]}, declination={self.declinations[0]}, cone={cone}, minloge={logE_threshold}")
                    lim = mean_limit(b, self.conflevel, False)
                    mrf = lim / s
                    mu_lds = get_mu_lds(b, self.significance, 0.5)
                    mdp = mu_lds / s
                    #print(f"Trying cone={cone}, minloge={logE_threshold} -> s={s}, b={b}, lim={lim}, mrf={mrf}, mu_lds={mu_lds}, mdp={mdp}")
                    return mrf
                
                for livetime in self.livetimes:  # Usa i parametri di livetime
                    for declination in self.declinations:  # Usa i parametri di declinazione
                        cone_start = 1.0
                        logE_threshold_start = 2.0

                        cone_range = [0.2, 5.0]
                        logE_threshold_range = [1., 3.5]

                        step_sizes = [0.05, 0.01]

                        # Ottimizzazione per il miglior cono e taglio logE
                        #cone, logE_threshold = pyroot_util.minimize(score, [1.0 , 0.5], [[0.2 , 5.0], [0., 8.0]], tollerance=1e-8)
                        
                        cone, logE_threshold = pyroot_util.minimize(
                                                score, 
                                                start_values=[cone_start, logE_threshold_start], 
                                                step_sizes=step_sizes,
                                                ranges=[cone_range, logE_threshold_range], 
                                                tollerance=1e-8
                                                )
                        
                        # Calcola i valori finali per il segnale e il rumore
                        #sig, bg, lim, mrf, mu_lds, mdp = self.get_sb_mod(flav, chan, livetime, declination, cone, minloge)
                        
                        s, b = self.get_sb_flav(flav, chan, livetime, declination, cone, logE_threshold)
                        
                        lim = mean_limit(b, self.conflevel, False)
                        mrf = lim / s
                        mu_lds = get_mu_lds(b, self.significance, 0.5)
                        mdp = mu_lds / s
                        # Aggiungi alla tabella
                        table.add(livetime).add(declination).add(cone ).add(logE_threshold).add(s).add(b).add(lim).add(mrf).add(mu_lds).add(mdp)

                        # Se l'output non è interattivo, scrivi su file
                        if self.output_style != "interact":
                            out_txt.write(f"{livetime} {declination} {s} {b} {lim} {mrf} {mu_lds} {mdp} {cone} {logE_threshold}\n")

        # Stampa la tabella
        print(table)

        # Chiudi il file se è stato aperto
        if self.output_style != "interact":
            out_txt.close()



    



if __name__ == "__main__":
    
    options  = aa.Options( __doc__.split('options:')[1] , sys.argv[1:] )


    print ( options )
    
    c = CalcRate(**aa.options_to_dict (options ))
    c.doit()
