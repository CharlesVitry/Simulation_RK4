from dataclasses import dataclass


@dataclass
class Population:
    TailleDeLaPopulationTotale: float
    Sains: float
    Infectes: float
    Retablis: float


@dataclass
class Parametre:
    JourDebut: int
    TauxRetablisement: float
    TauxInfection: float


@dataclass
class Scenario(Population):
    parametres: list[Parametre]
