# Hidraulikus hálózat részekre bontása

## Bevezetés

A `staci_split` program segítségével egy hidraulikus hálózat részhálózatokra bontható a következő szempontok szerint:

- Topológia: olyan részhálózatokat keresünk, melyeken belül maximális, ezeken kívül viszont minimális számú kapcsolat van.
- Érzékenység: ez esetben olyan csoportokat keresünk, melyeknek hasonló az össz-érzékenysége.
- Nyomásesés: ekkor a csoportok kumulált nyomásesése hasonló lesz.

## Quick start

A program futtatásához szükséges 

- a hidraulikus modellt tartalmazó `.spr` fájl, valamint
- a beállításokat tartalmazó `staci_split_settings.xml` fájl

A számítás maga a `staci_split` parancs kiadásával történik abban a könyvtárban, ahol a fenti két fájl is található.

A program futtatása után keletkező eredmények:

- `bog.dat`: az optimalizáló eljárás (genetikus algoritmus) kimenete (statisztika),
- `staci_split.log`: a staci_split üzenetei.
- `sensitivity_matrix_<...>.csv`: az érzékenységeket tartalmazó fájl.
- `membership.txt`: a csomópontok és ágak tagsága.

## A `staci_split_settings.xml` felépítése

Egy tipikus beállításokat tartalmazó fájl alább látható.

``` xml
<settings>
    <global_debug_level>2</global_debug_level>
    <n_comm>2</n_comm>
    <fname>10x10.spr</fname> 
    <obj_type>modularity</obj_type>
   <weight_type>sensitivity</weight_type>
    <weight_type_mod>friction_coeff</weight_type_mod>
    <logfilename> staci_split.log</logfilename>
    <popsize>50</popsize>
    <ngen>100</ngen>
    <pmut> 0.5</pmut>
    <pcross> 0.6</pcross>
</settings>
```

- `global_debug_level (0,1,2,3)`:  logfile részletessége, a növekvő érték egyre több információt jelent.
- `n_comm (int>0)`: a létrehozni kíván csoportok száma. Nagy érték megadása esetén előfordulhat, hogy a végül létrehozott csoport száma ennél kisebb, azaz keletkezne üres csoportok is.
- `fname (string)`: a staci projektfájl neve
- `obj_type (modularity|A-optimalty|D-optimality)`: az optimalizálás célfüggvénye. Csoportokra bontás esetén a `modularity` használandó
- `weight_tpye (friction_coeff|demand|diameter)`: a súlyozás lehetséges értékei:
	- `obj_type = modularity` esetén
		- `topology`:  topológia (klasszikus modularitás)
		- `dp`: nyomásesés alapján
		- `sensitivity`: érzékenységek, ez esetben a `weight_type_mod` mezőt `friction_coeff` (súrlódási tényező), `demand` (fogyasztás) vagy `diameter` (csőátmérő) értékre lehet állítani
	-  `obj_type = A-optimality | D-optimality` esetén
		- `friction_coeff | diameter | demand` lehet.
- `logfilename`: a program kimenete, a az információgazdagság a `global_debug_level` értékével állítható be.
- `popsize`: GA populáció mérete
- `ngen`: GA generációk száma
- `pmut`: GA mutációs ráta
- `pcross`: GA keresztezés valószínűsége