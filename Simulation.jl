  #ANNOTAZIONI:
#-Lettere greche: angoli in rad, lettere latine: angoli in gradi
#------------------------------------
using Stripeline
using AstroLib, Dates
using DelimitedFiles
using Printf
using Plots
using Healpix
using FITSIO
using DataFrames
using Rotations
using CSV

output_dir=raw"C:\Users\Leo\OneDrive - studenti.unimi.it\OneDrive - Università degli Studi di Milano\UNI\Tesi"
RPM = 1  #velocità motore zenitale
Ang_Lim = 10
Ang_Zenith = 20
start_day = DateTime( 2022, 1, 1, 0, 0)
telescope_motors(time_s) = (0.0, deg2rad(Ang_Zenith), timetorotang(time_s, RPM))
#Creo le struct dei vari corpi che analizzo
struct SUN end
struct Moon end
#DICTIONARY dei pianeti!
Planets=Dict("Mercury"=>1, "Venus"=>2, "Mars"=>4, "Jupiter"=>5, "Saturn"=>6, "Uranus"=>7, "Neptune"=>8, "Pluto"=>9)

function Horn_dir(H_ID::String, date, time_range_s) #Dove punta l'Horn H_ID a partire da data "date"
    db = InstrumentDB()
    H_o=db.focalplane[H_ID].orientation
    dirs, orientations = genpointings(
        telescope_motors,
        H_o,
        time_range_s,
        date)
    return dirs
end

function Horn_dir(H_o::Array, date, time_range_s)
    dirs, orientations = genpointings(
        telescope_motors,
        H_o,
        time_range_s,
        date)
    return dirs
end

function scalar_prod(v, w)
    aux = sum(v.*w)
end

function mod(v)
    aux = sqrt(sum(v.^2))
end

function ang_dist(α1, δ1, α2, δ2) #Distanza angolare
    θ = acos(sin(δ1)*sin(δ2) + cos(δ1)*cos(δ2)*cos(α1 - α2))
    return θ
end

function ang_dist(v, w) #Distanza angolare tra due vettori
    θ = acos(scalar_prod(v, w)/(mod(v)*mod(w)))
    return θ
end

function Pos(T::Type{SUN}, date)
    α, δ,_,_= AstroLib.sunpos(jdcnv(date), radians=true) #Da valori in rad
    return α, δ
end

function Pos(T::Type{Moon}, date)
    α, δ,_,_,_= AstroLib.moonpos(jdcnv(date), radians=true) #Da valori in rad
    return α, δ
end

function Pos(Pl_ID::Int64, date)
    ra, dec= AstroLib.planet_coords(jdcnv(date), Pl_ID) #Da valori in gradi
    α, δ = (deg2rad(ra), deg2rad(dec))
    return α, δ
end

function Pos(Pl_Name::String, date)
    Pl_ID = Planets[Pl_Name]
    ra, dec= AstroLib.planet_coords(jdcnv(date), Pl_ID) #Da valori in gradi
    α, δ = (deg2rad(ra), deg2rad(dec))
    return α, δ
end

function Dist(Obj, H_dir, date) #Restituisce posizione di un corpo celeste Obj e distanza
    α2, δ2 = Pos(Obj, date)     #dall'Horn con direzione H_dir in data date
    δ1, α1 = H_dir[1], H_dir[2] #dirs sono Dec e AR in rad.
    θ = ang_dist(α1, δ1, α2, δ2)
    return θ, α2, δ2
end

function Map_Plot_Horn(H_ID, start_day, sampling_time_s, lasting_time_s, map, save::Bool)
    filename=joinpath(output_dir,"$H_ID.map")
    time0=datetime2unix(start_day)
    time_range = 0:sampling_time_s:lasting_time_s
    db = InstrumentDB()
    H_o=db.focalplane[H_ID].orientation
    H_dirs=Horn_dir(H_o, start_day, time_range)
    for indx in 1:length(collect(time_range))
        pixel_index = ang2pix(map, π/2 - H_dirs[indx, 1], H_dirs[indx, 2])
        map[pixel_index] += 1 #aumento hit count del telescopio
    end
    if save==true
        saveToFITS(map, "!"*filename)
        return filename
    end
end

function Create_Map(Horns, lasting_time_s, name)
    #filename=joinpath(output_dir,"$name.map")
    map = Healpix.Map{Float64, RingOrder}(128)
    for indx in 1:length(Horns)
        Map_Plot_Horn(Horns[indx], DateTime( 2022, 1, 1, 0, 0), 0.1, lasting_time_s, map, false)
    end
    plot(map, mollweide, title="Map $name")
    saveToFITS(map, "!$name.map")
    return map, "$name.map"
end

function Angular_Ring(ax, θ, map, colour)  #plotta su mappa map anello di apertura θ con asse ax
    R1 = AngleAxis(θ, 1, 1, -(ax[1]+ax[2])/ax[3]) #trovo un qualunque ortogonale ad ax
    v0 = R1*ax   #e ruoto ax di θ attorno ad esso ottenendo v0
    for ϕ in range(0, 2*π, length=100)
        R2 = AngleAxis(ϕ, ax[1], ax[2], ax[3])
        v = R2*v0   #procedo a step ruotando v0 attorno all'asse
        δ, a = vec2ang(v[1], v[2], v[3])
        map[ang2pix(map, δ, a)] = colour #man mano coloro il pixel ottenuto
    end
end

function Map_Plot_Obj(Object, start_day, sampling_time_s, lasting_time_s, map, colour)
    time0=datetime2unix(start_day)
    time_range = 0:sampling_time_s:lasting_time_s
    time_samples = collect(time_range)
    for indx in 1:length(time_samples)
        time = time0 + indx*sampling_time_s
        α, δ= Pos(Object, unix2datetime(time))
        pixel_index_obj = ang2pix(map, π/2 - δ, α) #da declinazione (latitudine) a colat.
        map[pixel_index_obj] = colour #coloro la traccia dell'oggetto
    end
end

function Map_Plot(Object, H_ID, start_day, sampling_time_s, lasting_time_s)
    map = Healpix.Map{Float64, RingOrder}(128)
    Map_Plot_Horn(H_ID, start_day, sampling_time_s, lasting_time_s, map, false)
    max_hit=maximum(map)
    Map_Plot_Obj(Object, start_day, sampling_time_s, lasting_time_s, map, max_hit+max_hit/20)
    plot(map, mollweide, title="$Object - Horn $H_ID for $lasting_time_days days")
    savefig(joinpath(output_dir,"$Object-$H_ID-$lasting_time_days.pdf"))
end

#Overloading per poter passare la mappa dell'Horn già fatta!
function Map_Plot(Object, start_day, sampling_time_s, lasting_time_s, map_name)
    map = readMapFromFITS(map_name, 1, Float64)
    max_hit=maximum(map)
    Map_Plot_Obj(Object, start_day, sampling_time_s, lasting_time_s, map, max_hit+max_hit/20)
    lasting_time_days = lasting_time_s/86400
    plot(map, mollweide, title="$Object on $map_name for $lasting_time_days days")
    savefig(joinpath(output_dir,"$Object-$map_name-$lasting_time_days.pdf"))
end

#Overloading per stampare mappe senza oggetti sopra
function Map_Plot(map_name)
    map = readMapFromFITS(map_name, 1, Float64)
    plot(map, mollweide, title="$map_name")
    savefig("$map_name.pdf")
end

function Animation(Object, H_ID::String, start_day, lasting_time_weeks)
    map = Healpix.Map{Float64, RingOrder}(128)
    Map_Plot_Horn(H_ID, start_day, 11, lasting_time_weeks*7*86400, map, false)
    max_hit=maximum(map)
    sampling_time_s = 3600
    anim = @animate for week in 0:lasting_time_weeks
        start_time_s = datetime2unix(start_day) + week*7*86400
        one_week_s = 86400*7 - sampling_time_s #tolgo l'ultimo sample perchè verrà preso dalla settimana dopo a t=0
        Map_Plot_Obj(Object, unix2datetime(start_time_s), sampling_time_s, one_week_s, map, max_hit+max_hit/20)
        plot(map, mollweide,title="$Object - Horn $H_ID week $week")
    end
    gif(anim, joinpath(output_dir,"$Object-$H_ID-animation.gif"), fps = 10) #velocità
end

#Creo un overloading di Animation, in modo da potergli anche passare la mappa (salvata in FITS)
#con la posizione dell'Horn già plottata. Mi permette di creare una volta sola la
#mappa statica dell'Horn molto più dettagliata e farci sopra animation per
#corpi diversi in modo molto più veloce. Non ha più bisogno di H_ID in imput.
function Animation(Object, start_day::DateTime, lasting_time_weeks, map_name::String)
    map = readMapFromFITS(map_name, 1, Float64)
    max_hit=maximum(map)
    sampling_time_s = 3600
    anim = @animate for week in 0:lasting_time_weeks
        start_time_s = datetime2unix(start_day) + week*7*86400
        one_week_s = 86400*7 - sampling_time_s #tolgo l'ultimo sample perchè verrà preso dalla settimana dopo a t=0
        Map_Plot_Obj(Object, unix2datetime(start_time_s), sampling_time_s, one_week_s, map, max_hit+max_hit/20)
        α, δ = Pos(Object, unix2datetime(start_time_s + one_week_s))
        ax = ang2vec(π/2 - δ, α)
        Angular_Ring(ax, deg2rad(10), map, max_hit+max_hit/20)
        plot(map, mollweide, title="$Object on $map_name week $week")
        Angular_Ring(ax, deg2rad(10), map, 0)
    end
    gif(anim, joinpath(output_dir,"$Object-$map_name-animation.gif"), fps = 10) #velocità
end

function Animation_Sun(Object, start_day::DateTime, lasting_time_weeks, map_name::String)
    map = readMapFromFITS(map_name, 1, Float64)
    max_hit=maximum(map)
    sampling_time_s = 3600
    anim = @animate for week in 0:lasting_time_weeks
        start_time_s = datetime2unix(start_day) + week*7*86400
        one_week_s = 86400*7 - sampling_time_s #tolgo l'ultimo sample perchè verrà preso dalla settimana dopo a t=0
        Map_Plot_Obj(Object, unix2datetime(start_time_s), sampling_time_s, one_week_s, map, max_hit+max_hit/20)
        Map_Plot_Obj(SUN, unix2datetime(start_time_s), sampling_time_s, one_week_s, map, max_hit*2/3)
        α, δ = Pos(SUN, unix2datetime(start_time_s + one_week_s))
        ax = ang2vec(π/2 - δ, α)
        Angular_Ring(ax, deg2rad(10), map, max_hit*2/3)
        plot(map, mollweide, title="$Object + Sun on $map_name week $week")
        map = readMapFromFITS(map_name, 1, Float64)
    end
    gif(anim, joinpath(output_dir,"$Object-Sun-$map_name-animation.gif"), fps = 10) #velocità
end

function Animation_Horns(start_day::DateTime, lasting_time_s)
    map1 = Healpix.Map{Float64, RingOrder}(128)
    start_time_s = datetime2unix(start_day)
    time_range = 0:0.5:lasting_time_s
    time_samples = collect(time_range)
    G6_dirs=Horn_dir("G6", start_day, time_samples)
    I0_dirs=Horn_dir("I0", start_day, time_samples)
    V3_dirs=Horn_dir("V3", start_day, time_samples)
    anim = @animate for indx in 1:length(time_samples)
        pixel_index = ang2pix(map1, π/2 - G6_dirs[indx, 1], G6_dirs[indx, 2])
        map1[pixel_index] += 1
        pixel_index = ang2pix(map1, π/2 - I0_dirs[indx, 1], I0_dirs[indx, 2])
        map1[pixel_index] += 1
        pixel_index = ang2pix(map1, π/2 - V3_dirs[indx, 1], V3_dirs[indx, 2])
        map1[pixel_index] += 1
        plot(map1, mollweide)
    end
    gif(anim, joinpath(output_dir,"Telescope_animation.gif"), fps=20)
end

#time0 è in formato UNIX!!!
function OneH_Simulation(Object, H_o, time0, sampling_time_s, lasting_time_s, f)
    time_range = 0:sampling_time_s:lasting_time_s
    time_samples = collect(time_range)
    H_dirs=Horn_dir(H_o, start_day, time_samples)
    for indx in 1:length(time_samples)
        time = time0 + indx*sampling_time_s
        date=unix2datetime(time)
        θ, α, δ=Dist(Object, H_dirs[indx, :], date)
        if θ < deg2rad(20)
            @printf(f, "%s,%.5e,%.5e,%.5e\n", date, θ, α, δ)
        end
    end
end

function Simulation(Object, H_ID::String, start_day)
    db = InstrumentDB()
    H_o=db.focalplane[H_ID].orientation
    sampling_time_s=0.1 #va messo 0.1
    lasting_time_s=86400
    time0=datetime2unix(start_day)
    filename=joinpath(output_dir,"$Object.$H_ID.txt")
    open(filename, "w") do f
        writedlm(f, ["date" "Ang. Dist. $Object-Horn $H_ID [rad]" "$Object's RA [rad]" "$Object's Dec [rad]"], ',')
        for day in 1:730 #2 anni = 730giorni
            OneH_Simulation(Object, H_o, time0+86400*day, sampling_time_s, lasting_time_s, f)
        end
    end
    filename
end

function Simulation(Object, H_o::Array, start_day)
    sampling_time_s=0.1 #va messo 0.1
    lasting_time_s=86400
    time0=datetime2unix(start_day)
    filename=joinpath(output_dir,"$Object.$H_ID.txt")
    open(filename, "w") do f
        writedlm(f, ["date" "Ang. Dist. $Object-Horn $H_ID [rad]" "$Object's RA [rad]" "$Object's Dec [rad]"], ',')
        for day in 1:730 #2 anni = 730giorni
            OneH_Simulation(Object, H_o, time0+86400*day, sampling_time_s, lasting_time_s, f)
        end
    end
    filename
end

function File2dates(filename, ang_lim)
    dat = readdlm(filename,  ',', skipstart=1)
    dates = dat[dat[:,2] .< deg2rad(ang_lim),1]
    return dates
end

function Dates_Union(dates_matr)
    dates = dates_matr[1]
    for indx in 2:length(dates_matr)
        dates_aux = dates_matr[indx]
        dates = union(dates, dates_aux) #faccio l'unione delle date
    end
    sort!(dates)
    return dates
end

function Dates_Intersection(dates_matr)
    dates = dates_matr[1]
    for indx in 2:length(dates_matr)
        dates_aux = dates_matr[indx]
        dates = intersect(dates, dates_aux) #faccio l'intersezione delle date
    end
    sort!(dates)
    return dates
end

function Time_Count(dates)
    rot_T_s = 60/RPM
    len = length(dates[:,1])
    time_s = 0
    i=1
    while i < len
        t0 = datetime2unix(DateTime(dates[i]))
        while datetime2unix(DateTime(dates[i]))-t0 < rot_T_s && i < len
            i += 1
        end
        time_s += rot_T_s
    end
    return time_s
end

function Planet_Result(files, ang_lims, sampling_time_s, output_name)
    filename = joinpath(output_dir,output_name*".txt")
    start_day_u=datetime2unix(start_day)
    times = []
    dates_matr = []
    for (indx, file) in enumerate(files)  #estraggo le date dai vari file
        dates = File2dates(file, ang_lims[indx])
        push!(dates_matr, dates)
    end
    dates = Dates_Union(dates_matr)    #faccio l'unione
    L_0 = 0
    for week in 1:104   #calcolo i giorni di osservazione in una settimana
        w_end_u = start_day_u + week*7*86400
        L = length(dates[dates .< string(unix2datetime(w_end_u))])
        t = L - L_0
        push!(times, t*sampling_time_s/3600)  #converto in ore
        L_0 = L
    end
    open(filename, "w") do f
        writedlm(f, times)
    end
    return times
end

function Planets_Plot(Planets)
    filename = joinpath(output_dir, "Planets.pdf")
    files = Planets.*"_results.txt"
    plot(1:104, readdlm(files[1]), label = Planets[1], title="Planets Observed ($Ang_Lim ° threshold)", xlabel="Weeks since 01/01/2022", ylabel="Observation time [h/week]")
    for indx in 2:length(files)
        times = readdlm(files[indx])
        plot!(1:104, times, label = Planets[indx])
    end
    savefig(filename)
end

function Get_Results(files, ang_lims, output_name)
    filename = joinpath(output_dir,output_name*".txt")
    times = []
    dates_matr = []
    for (indx, file) in enumerate(files) #faccio il conto dei singoli Horn
        dates = File2dates(file, ang_lims[indx])
        push!(dates_matr, dates)
        t_s = Time_Count(dates)
        push!(times, t_s/86400)
    end
    push!(files, "total")  #faccio il conto totale unendo le info
    push!(ang_lims, Ang_Lim)
    t_s = Time_Count(Dates_Union(dates_matr))
    push!(times, t_s/86400)
    percs = times./730*100
    df = DataFrame(file=files, ang_lim=ang_lims, time_days=times, percentage=percs)
    open(filename, "w") do f
    show(f, MIME("text/latex"), df, eltypes=false)
    end
end

#-------------------------   "MAIN"
#"!" * joinpath(percorso, "$H_ID.map")

files = ["SUN.G6.txt", "SUN.I0.txt", "SUN.V3.txt", "Moon.G6.txt", "Moon.I0.txt", "Moon.V3.txt", "total"]
ang_lims = [Ang_Lim, Ang_Lim, Ang_Lim, Ang_Lim, Ang_Lim, Ang_Lim, Ang_Lim]
times = [58.8361, 50.8889, 42.7972, 58.8965, 51.8431, 44.5215, 135.458]
percs = [8.05974, 6.97108, 5.86263, 8.06802, 7.10179, 6.09884, 18.5559]
df = DataFrame(file=files, ang_lim=ang_lims, time_days=times, percentage=percs)
open("total_latex.txt", "w") do f
    show(f, MIME("text/latex"), df, eltypes=false)
end

#SIMULAZIONE
Object="Neptune"
Simulation(Object, "G6", start_day)
Simulation(Object, "I0", start_day)
Simulation(Object, "V3", start_day)
Simulation(Object, "B5", start_day)

H_files=Array(["Neptune.G6.txt", "Neptune.I0.txt", "Neptune.V3.txt"])
angs=repeat(Ang_Lim:Ang_Lim, outer=length(H_files))
Get_Results(H_files, angs, "Moon_results")

H_files=Array(["SUN.G6.txt"])
angs=Array([Ang_Lim])
Get_Results(H_files, angs, "prova")

Planet_Result(H_files, angs, 0.1, Object*"_results")
Planets_Plot(Array([ "Mercury", "Venus", "Mars", "Jupiter", "Uranus", "Neptune"]))

start_date = DateTime(datas[2, 1])
plot([(DateTime(x) - start_date).value*1e-3/86400 for x in datas[2:length(datas[:,1]), 1]], ones(length(datas[:,1])-1))
#rm("$Object.$Horn.txt")

filename = joinpath(output_dir, "Sun_Moon.pdf")
plot(1:104, Planet_Result(["SUN.G6.txt", "SUN.I0.txt", "SUN.V3.txt"], [10, 10, 10], 0.1, "Sun_res"), label = "Sun", title="Sun/Moon interference with observations ($Ang_Lim ° threshold)", xlabel="Weeks since 01/01/2022", ylabel="Time discarded [h/week]")
plot!(1:104, Planet_Result(["Moon.G6.txt", "Moon.I0.txt", "Moon.V3.txt"], [10, 10, 10], 0.1, "Moon_res"), label = "Moon")
savefig(filename)


#ANIMAZIONI
#Creo la mappa per un Horn, 1 volta sola!! (su 20gg è già abbastanza definita)
lasting_time_days = 1
lasting_time_s=3600
Horns = ["I0"]
m, _ = Create_Map(Horns, lasting_time_s, "I0-1h")
plot(m, mollweide)

#Ci faccio sopra l' animazione del corpo celeste che voglio
Animation_Sun("Mercury", start_day, 104, "G6+I0+V3.map")

#ANIMAZIONE TELESCOPIO
Animation_Horns(start_day, 120)

#PLOT STATICI
Map_Plot("Saturn", start_day, 3600, 86400*365*2, "G6+I0+V3.map")
savefig(joinpath(output_dir,"$Object-$Horn-$lasting_time_days.pdf"))

Map_Plot("I0-1h.map")

map1 = Healpix.Map{Float64, RingOrder}(128)
Angular_Ring(ang2vec(π/2, π/2), deg2rad(10), map1, 4)
plot(map1, mollweide)

#Calcolo horn "fittizio"
db = InstrumentDB()
#i 3 horn che mi interessano
I0 = db.focalplane["I0"].orientation
B0 = db.focalplane["B0"].orientation
O0 = db.focalplane["O0"].orientation
#mi metto nel SR del boresight, lo zenith è 20° indietro nella direzione dell x (asse di simmetria del telesc)
zenith = [sin(deg2rad(-Ang_Zenith)), 0, cos(deg2rad(-Ang_Zenith))]
#calcolo δ: la differenza di distanza tra O0 e lo zenith (con B0 è uguale per simmetria) e I0 e lo zenith
δ = ang_dist(zenith, O0) - ang_dist(zenith, I0)
#creo delle coordinate di un fake_horn nella stessa posizione di I0 ma ruotato di δ
#per avere lo stesso angolo con lo zenith di B0 e O0.
f_h = [sin(δ), 0, cos(δ)]
#posso passare f_h al posto di H_ID grazie ad un overloading di H_dir.

#I1 I5
#O1 B5
