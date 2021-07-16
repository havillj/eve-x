import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from lxml import etree
from matplotlib.backends.backend_pdf import PdfPages
import io
import base64
from dna_features_viewer import GraphicFeature, GraphicRecord
from dash.dependencies import Input, Output, State
import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.use('Agg')
import math
import dash_bootstrap_components as dbc


specVirus = {('Philippines', 15) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Merida virus', 'Wuhan Mosquito Virus 6', 'Zika virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'West Nile virus'] ,
('Angola', 9) : ['Lassa mammarenavirus', 'Liao ning virus', 'Alfalfa mosaic virus', 'Ohlsdorf ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Zika virus', 'Radi vesiculovirus', 'Phasi Charoen-like phasivirus'] ,
('Brazil', 2) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Alfalfa mosaic virus', 'Phasi Charoen-like phasivirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Avian orthoavulavirus 1', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Liao ning virus'] ,
('Philippines', 6) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Prunus necrotic ringspot virus', 'Pestivirus A', 'Australian Anopheles totivirus', 'West Nile virus', 'Aedes anphevirus', 'Riverside virus 1', 'Cell fusing agent virus', 'Menghai flavivirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Gabon', 6) : ['Alfalfa mosaic virus', 'Pestivirus A', 'Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Phasi Charoen-like phasivirus', 'Xishuangbanna aedes flavivirus', 'Ohlsdorf ohlsrhavirus', 'Ochlerotatus scapularis flavivirus'] ,
('Angola', 1) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Avian orthoavulavirus 1', 'Aedes aegypti anphevirus', 'Tongilchon ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Zika virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Ochlerotatus scapularis flavivirus'] ,
('Angola', 16) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Aedes aegypti anphevirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Zika virus', 'Riverside virus 1', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Culex Flavi-like virus', 'Wuhan Mosquito Virus 6'] ,
('Argentina', 2) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Fitzroy Crossing toti-like virus 1', 'Culex ohlsrhavirus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Alfalfa mosaic virus', 'Culex flavivirus'] ,
('French-Polynesia', 9) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Tongilchon ohlsrhavirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus'] ,
('South-Africa', 14) : ['Mercadeo virus', 'Prunus necrotic ringspot virus', 'Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Phasi Charoen-like phasivirus', 'Zika virus'] ,
('Australia', 10) : ['Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Liao ning virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Tongilchon ohlsrhavirus', 'Zika virus', 'Xishuangbanna aedes flavivirus'] ,
('Mexico', 24) : ['Lassa mammarenavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Ochlerotatus scapularis flavivirus', 'Merida virus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Avian orthoavulavirus 1', 'Pestivirus A'] ,
('USA', 20) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Prunus necrotic ringspot virus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Pestivirus A'] ,
('Brazil', 19) : ['Merida virus', 'Liao ning virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Lassa mammarenavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Alfalfa mosaic virus', 'Culex ohlsrhavirus', 'Culex flavivirus', 'Zika virus'] ,
('Gabon', 3) : ['Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Merida virus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Avian orthoavulavirus 1'] ,
('Gabon', 16) : ['Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Merida virus', 'Wuhan Mosquito Virus 6', 'Radi vesiculovirus', 'Culex ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Mercadeo virus', 'Tongilchon ohlsrhavirus'] ,
('Mexico', 18) : ['Lassa mammarenavirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Mercadeo virus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Thailand', 14) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Pestivirus A', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Zika virus', 'Ochlerotatus scapularis flavivirus'] ,
('Angola', 10) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Liao ning virus', 'Zika virus', 'Phasi Charoen-like phasivirus', 'Ohlsdorf ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Avian orthoavulavirus 1'] ,
('Australia', 1) : ['Wuhan Mosquito Virus 6', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Alfalfa mosaic virus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Zika virus', 'Pestivirus A'] ,
('USA', 9) : ['Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Tongilchon ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Mercadeo virus', 'Pestivirus A'] ,
('Australia', 20) : ['Australian Anopheles totivirus', 'Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Avian orthoavulavirus 1', 'Zika virus', 'Tick-borne encephalitis virus', 'Pestivirus A'] ,
('Vietnam', 2) : ['Mercadeo virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Zika virus'] ,
('Thailand', 7) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Ochlerotatus scapularis flavivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Zika virus', 'Alfalfa mosaic virus', 'Menghai flavivirus', 'Australian Anopheles totivirus'] ,
('Philippines', 22) : ['Liao ning virus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Xishuangbanna aedes flavivirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Mercadeo virus', 'Pestivirus A', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Zika virus', 'Ochlerotatus scapularis flavivirus', 'Prunus necrotic ringspot virus'] ,
('Gabon', 10) : ['Pestivirus A', 'Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'Phasi Charoen-like phasivirus'] ,
('Mexico', 16) : ['Lassa mammarenavirus', 'Parramatta River virus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Phasi Charoen-like phasivirus', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Zika virus'] ,
('Australia', 12) : ['Lassa mammarenavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Avian orthoavulavirus 1', 'Ochlerotatus scapularis flavivirus', 'Zika virus', 'Pestivirus A'] ,
('Gabon', 19) : ['Menghai flavivirus', 'Lassa mammarenavirus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Riverside virus 1', 'Culiseta flavivirus', 'Pestivirus A', 'West Nile virus', 'Zika virus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6'] ,
('South-Africa', 12) : ['Lassa mammarenavirus', 'Culex Flavi-like virus', 'Wuhan Mosquito Virus 6', 'Alfalfa mosaic virus', 'Prunus necrotic ringspot virus', 'Ohlsdorf ohlsrhavirus', 'Radi vesiculovirus', 'Australian Anopheles totivirus', 'Zika virus', 'Phasi Charoen-like phasivirus', 'Aedes anphevirus', 'Tongilchon ohlsrhavirus'] ,
('USA', 7) : ['Prunus necrotic ringspot virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus', 'West Nile virus', 'Zika virus'] ,
('Vietnam', 9) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Fitzroy Crossing toti-like virus 1', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Pestivirus A', 'Menghai flavivirus', 'Zika virus'] ,
('Philippines', 23) : ['Alfalfa mosaic virus', 'Culiseta flavivirus', 'Menghai flavivirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Avian orthoavulavirus 1', 'West Nile virus', 'Zika virus', 'Lassa mammarenavirus'] ,
('Brazil', 24) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Culex ohlsrhavirus'] ,
('Mexico', 20) : ['Lassa mammarenavirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Phasi Charoen-like phasivirus', 'Cell fusing agent virus', 'Ochlerotatus caspius flavivirus', 'Ochlerotatus scapularis flavivirus'] ,
('French-Polynesia', 18) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Xishuangbanna aedes flavivirus', 'Human orthopneumovirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Mexico', 11) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Ohlsdorf ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Avian orthoavulavirus 1'] ,
('Thailand', 17) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Zika virus', 'Alfalfa mosaic virus'] ,
('Gabon', 5) : ['Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'West Nile virus', 'Grenada mosquito rhabdovirus 1', 'Mercadeo virus'] ,
('Mexico', 5) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus', 'Avian orthoavulavirus 1'] ,
('Angola', 12) : ['Lassa mammarenavirus', 'Liao ning virus', 'Ohlsdorf ohlsrhavirus', 'Alfalfa mosaic virus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Zika virus', 'Wuhan Mosquito Virus 6', 'Culiseta flavivirus', 'Gata virus'] ,
('Mexico', 2) : ['Lassa mammarenavirus', 'Liao ning virus', 'Mercadeo virus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Pestivirus A', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Menghai flavivirus', 'Zika virus', 'Avian orthoavulavirus 1', 'Riverside virus 1'] ,
('Philippines', 21) : ['Mercadeo virus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Pestivirus A', 'Lassa mammarenavirus', 'West Nile virus', 'Zika virus'] ,
('Argentina', 3) : ['Mercadeo virus', 'Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Menghai flavivirus', 'Zika virus', 'Xishuangbanna aedes flavivirus'] ,
('USA', 10) : ['Lassa mammarenavirus', 'Merida virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Phasi Charoen-like phasivirus', 'Ochlerotatus caspius flavivirus', 'Xishuangbanna aedes flavivirus', 'Cell fusing agent virus', 'Liao ning virus', 'Zika virus'] ,
('Gabon', 11) : ['Australian Anopheles totivirus', 'Lassa mammarenavirus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Phasi Charoen-like phasivirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Mac Peak virus', 'Pestivirus A', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('USA', 1) : ['Lassa mammarenavirus', 'Parramatta River virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Pestivirus A', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Ohlsdorf ohlsrhavirus', 'Mansonia flavivirus', 'Grenada mosquito rhabdovirus 1'] ,
('USA', 4) : ['Merida virus', 'Xishuangbanna aedes flavivirus', 'Kamiti River virus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Lassa mammarenavirus', 'Zika virus'] ,
('Gabon', 18) : ['Cell fusing agent virus', 'Lassa mammarenavirus', 'Tongilchon ohlsrhavirus', 'Culiseta flavivirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Merida virus', 'Ohlsdorf ohlsrhavirus', 'Culex ohlsrhavirus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Argentina', 8) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Prunus necrotic ringspot virus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Australian Anopheles totivirus', 'Riverside virus 1', 'Culex ohlsrhavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Avian orthoavulavirus 1', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'West Nile virus', 'Menghai flavivirus'] ,
('Vietnam', 16) : ['Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Mercadeo virus', 'Aedes anphevirus', 'Liao ning virus', 'Menghai flavivirus', 'Lassa mammarenavirus', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Culiseta flavivirus'] ,
('USA', 23) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Pestivirus A', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Philippines', 12) : ['Mercadeo virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Lassa mammarenavirus', 'Zika virus', 'Pestivirus A'] ,
('French-Polynesia', 15) : ['Cell fusing agent virus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Culiseta flavivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Thailand', 15) : ['Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Lassa mammarenavirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Ochlerotatus caspius flavivirus', 'Pestivirus A', 'Quang Binh virus'] ,
('French-Polynesia', 11) : ['Aedes anphevirus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Riverside virus 1', 'Liao ning virus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Alfalfa mosaic virus'] ,
('Brazil', 15) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Merida virus', 'West Nile virus', 'Cell fusing agent virus', 'Phasi Charoen-like phasivirus', 'Xishuangbanna aedes flavivirus', 'Mercadeo virus', 'Ochlerotatus caspius flavivirus', 'Australian Anopheles totivirus', 'Zika virus'] ,
('French-Polynesia', 14) : ['Lassa mammarenavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Fitzroy Crossing toti-like virus 1', 'Zika virus', 'Australian Anopheles totivirus', 'Xishuangbanna aedes flavivirus', 'Avian orthoavulavirus 1', 'Alfalfa mosaic virus'] ,
('French-Polynesia', 6) : ['Radi vesiculovirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Lassa mammarenavirus', 'Liao ning virus', 'Riverside virus 1', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Fitzroy Crossing toti-like virus 1', 'Ochlerotatus scapularis flavivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Culex ohlsrhavirus'] ,
('Gabon', 24) : ['Lassa mammarenavirus', 'Tongilchon ohlsrhavirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Culex mosquito virus 2', 'Ohlsdorf ohlsrhavirus', 'Australian Anopheles totivirus', 'Avian orthoavulavirus 1', 'Ochlerotatus scapularis flavivirus', 'Zika virus'] ,
('Gabon', 12) : ['Lassa mammarenavirus', 'Culiseta flavivirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Phasi Charoen-like phasivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Tongilchon ohlsrhavirus', 'Zika virus'] ,
('French-Polynesia', 10) : ['Lassa mammarenavirus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Aedes anphevirus', 'Aedes aegypti anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Avian orthoavulavirus 1', 'Riverside virus 1', 'Culex flavivirus', 'Ochlerotatus scapularis flavivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Liao ning virus', 'Alfalfa mosaic virus'] ,
('Vietnam', 3) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Liao ning virus', 'Ohlsdorf ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Kamiti River virus', 'Pestivirus A', 'Zika virus'] ,
('Philippines', 8) : ['Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Lassa mammarenavirus', 'Merida virus', 'Liao ning virus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Pestivirus A', 'Zika virus'] ,
('Vietnam', 17) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Mercadeo virus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Australia', 18) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Zika virus', 'Ohlsdorf ohlsrhavirus'] ,
('Thailand', 1) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Fitzroy Crossing toti-like virus 1', 'Pestivirus A', 'Zika virus'] ,
('Philippines', 2) : ['Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Pestivirus A', 'Cell fusing agent virus', 'Tongilchon ohlsrhavirus', 'Ochlerotatus caspius flavivirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Vietnam', 1) : ['Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Fitzroy Crossing toti-like virus 1', 'Liao ning virus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Menghai flavivirus'] ,
('Angola', 17) : ['Zika virus', 'Liao ning virus', 'Ochlerotatus caspius flavivirus', 'Radi vesiculovirus', 'Riverside virus 1', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Aedes aegypti anphevirus', 'Lassa mammarenavirus', 'Ohlsdorf ohlsrhavirus'] ,
('Thailand', 11) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Zika virus'] ,
('Australia', 17) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Tongilchon ohlsrhavirus', 'West Nile virus', 'Zika virus'] ,
('Australia', 16) : ['Liao ning virus', 'Alfalfa mosaic virus', 'Aedes anphevirus', 'Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Ochlerotatus scapularis flavivirus', 'Tongilchon ohlsrhavirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Zika virus', 'Wuhan Mosquito Virus 6', 'Drosophila melanogaster sigmavirus', 'Xishuangbanna aedes flavivirus'] ,
('Angola', 14) : ['Lassa mammarenavirus', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Phasi Charoen-like phasivirus', 'Culex ohlsrhavirus', 'Zika virus', 'Tongilchon ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Fitzroy Crossing toti-like virus 1', 'Australian Anopheles totivirus'] ,
('Australia', 11) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Alfalfa mosaic virus', 'Liao ning virus', 'Aedes anphevirus', 'Tongilchon ohlsrhavirus', 'Fitzroy Crossing toti-like virus 1', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Malpais Spring vesiculovirus'] ,
('USA', 2) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Xishuangbanna aedes flavivirus', 'Liao ning virus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Gabon', 17) : ['Lassa mammarenavirus', 'Dengue virus', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Xishuangbanna aedes flavivirus', 'Pestivirus A', 'Avian orthoavulavirus 1'] ,
('USA', 24) : ['Lassa mammarenavirus', 'West Nile virus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Aedes aegypti anphevirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Pestivirus A', 'Radi vesiculovirus', 'Ochlerotatus scapularis flavivirus', 'Zika virus'] ,
('Vietnam', 8) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Liao ning virus', 'Ohlsdorf ohlsrhavirus', 'Zika virus', 'Australian Anopheles totivirus', 'Ochlerotatus scapularis flavivirus', 'Prunus necrotic ringspot virus'] ,
('Brazil', 23) : ['Lassa mammarenavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Zika virus', 'Pestivirus A'] ,
('USA', 16) : ['Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Ochlerotatus scapularis flavivirus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Pestivirus A', 'Zika virus'] ,
('French-Polynesia', 12) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Ohlsdorf ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes aegypti anphevirus', 'Liao ning virus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Riverside virus 1', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Pestivirus A', 'Zika virus'] ,
('USA', 19) : ['Lassa mammarenavirus', 'Liao ning virus', 'Merida virus', 'Ochlerotatus caspius flavivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'West Nile virus', 'Tongilchon ohlsrhavirus', 'Culex ohlsrhavirus', 'Avian orthoavulavirus 1', 'Zika virus'] ,
('Brazil', 13) : ['Liao ning virus', 'Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Aedes aegypti anphevirus', 'Cell fusing agent virus', 'Phasi Charoen-like phasivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ochlerotatus scapularis flavivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Aedes anphevirus'] ,
('USA', 5) : ['Lassa mammarenavirus', 'Merida virus', 'Aedes anphevirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'West Nile virus', 'Japanese encephalitis virus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus'] ,
('Mexico', 1) : ['Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Aedes flavivirus', 'Pestivirus A'] ,
('French-Polynesia', 2) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Liao ning virus', 'Tongilchon ohlsrhavirus', 'Fitzroy Crossing toti-like virus 1', 'Wuhan Mosquito Virus 6', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Xishuangbanna aedes flavivirus'] ,
('Gabon', 9) : ['Lassa mammarenavirus', 'Culiseta flavivirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Tongilchon ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Zika virus'] ,
('Philippines', 4) : ['Lassa mammarenavirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Mercadeo virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'Pestivirus A', 'Zika virus'] ,
('USA', 15) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Philippines', 5) : ['Liao ning virus', 'Merida virus', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Pestivirus A', 'West Nile virus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Philippines', 20) : ['Wuhan Mosquito Virus 6', 'Culiseta flavivirus', 'Lassa mammarenavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Zika virus'] ,
('Angola', 5) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Ochlerotatus scapularis flavivirus', 'Aedes anphevirus', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Zika virus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes aegypti anphevirus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'Avian orthoavulavirus 1'] ,
('French-Polynesia', 20) : ['Lassa mammarenavirus', 'Ohlsdorf ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Fitzroy Crossing toti-like virus 1', 'Culex tritaeniorhynchus rhabdovirus', 'La Tina virus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Culex ohlsrhavirus'] ,
('Philippines', 10) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Pestivirus A', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Zika virus', 'Avian orthoavulavirus 1', 'West Nile virus'] ,
('Vietnam', 19) : ['Lassa mammarenavirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Zika virus'] ,
('French-Polynesia', 5) : ['Wuhan Mosquito Virus 6', 'Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Aedes anphevirus', 'Tongilchon ohlsrhavirus', 'Liao ning virus', 'Australian Anopheles totivirus', 'Radi vesiculovirus', 'Cell fusing agent virus', 'Avian orthoavulavirus 1', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Brazil', 1) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Parramatta River virus'] ,
('Philippines', 1) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Riverside virus 1', 'Cell fusing agent virus', 'Mercadeo virus', 'Avian orthoavulavirus 1', 'Pestivirus A', 'Zika virus', 'Xingshan alphanemrhavirus', 'Ochlerotatus caspius flavivirus'] ,
('Philippines', 14) : ['Lassa mammarenavirus', 'Liao ning virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Mercadeo virus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Avian orthoavulavirus 1', 'Zika virus'] ,
('Angola', 8) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Liao ning virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Aedes aegypti anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Zika virus', 'Mercadeo virus', 'Riverside virus 1'] ,
('Brazil', 6) : ['Aedes anphevirus', 'Tongilchon ohlsrhavirus', 'Lassa mammarenavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Alfalfa mosaic virus', 'Pestivirus A', 'Australian Anopheles totivirus', 'Zika virus', 'Avian orthoavulavirus 1', 'Serbia mononega-like virus 1'] ,
('Mexico', 10) : ['Liao ning virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus', 'La Tina virus', 'Culiseta flavivirus', 'Lassa mammarenavirus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Mexico', 22) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Aedes aegypti anphevirus', 'Tongilchon ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'West Nile virus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Avian orthoavulavirus 1'] ,
('Philippines', 16) : ['American plum line pattern virus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'West Nile virus', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Tongilchon ohlsrhavirus', 'Pestivirus A', 'Mercadeo virus', 'Zika virus'] ,
('Brazil', 7) : ['Lassa mammarenavirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Mercadeo virus', 'Liao ning virus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Brazil', 4) : ['Ochlerotatus scapularis flavivirus', 'Culex ohlsrhavirus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Alfalfa mosaic virus', 'Ohlsdorf ohlsrhavirus', 'Zika virus', 'Australian Anopheles totivirus'] ,
('Vietnam', 15) : ['Merida virus', 'Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ohlsdorf ohlsrhavirus', 'Australian Anopheles totivirus', 'Pestivirus A', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Mexico', 8) : ['Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Liao ning virus', 'Ohlsdorf ohlsrhavirus', 'Avian orthoavulavirus 1', 'Wuhan Mosquito Virus 6', 'Ochlerotatus scapularis flavivirus', 'Phasi Charoen-like phasivirus', 'Zika virus', 'Australian Anopheles totivirus', 'Pestivirus A'] ,
('Mexico', 3) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Cell fusing agent virus', 'La Tina virus', 'Pestivirus A', 'Avian orthoavulavirus 1', 'Tongilchon ohlsrhavirus'] ,
('USA', 18) : ['Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'West Nile virus', 'Cell fusing agent virus', 'Lassa mammarenavirus', 'Merida virus', 'Liao ning virus', 'Phasi Charoen-like phasivirus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus', 'Aedes aegypti anphevirus', 'Avian orthoavulavirus 1'] ,
('Vietnam', 18) : ['Mercadeo virus', 'Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Zika virus', 'Culex tritaeniorhynchus rhabdovirus', 'Fitzroy Crossing toti-like virus 1', 'West Nile virus', 'Xishuangbanna aedes flavivirus', 'Menghai flavivirus'] ,
('Vietnam', 4) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Zika virus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus'] ,
('Thailand', 13) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Radi vesiculovirus', 'Zika virus', 'Avian orthoavulavirus 1', 'Xishuangbanna aedes flavivirus'] ,
('Thailand', 5) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Argentina', 11) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Merida virus', 'Liao ning virus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Pestivirus A', 'Aedes flavivirus'] ,
('Australia', 7) : ['Wuhan Mosquito Virus 6', 'Alfalfa mosaic virus', 'Australian Anopheles totivirus', 'Tongilchon ohlsrhavirus', 'Lassa mammarenavirus', 'Cell fusing agent virus', 'Liao ning virus', 'Aedes anphevirus', 'Menghai flavivirus', 'Zika virus', 'Parramatta River virus', 'Avian orthoavulavirus 1', 'Aedes aegypti anphevirus'] ,
('South-Africa', 9) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Zika virus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Australian Anopheles totivirus', 'Riverside virus 1', 'Radi vesiculovirus', 'Aedes anphevirus', 'Ochlerotatus scapularis flavivirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Cell fusing agent virus'] ,
('USA', 12) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Phasi Charoen-like phasivirus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Gabon', 7) : ['Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Merida virus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Ochlerotatus scapularis flavivirus', 'Avian orthoavulavirus 1', 'Culex tritaeniorhynchus rhabdovirus', 'Culiseta flavivirus', 'Radi vesiculovirus', 'Phasi Charoen-like phasivirus', 'Pestivirus A'] ,
('USA', 22) : ['Lassa mammarenavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Mercadeo virus', 'Xishuangbanna aedes flavivirus'] ,
('Gabon', 14) : ['Menghai flavivirus', 'Mercadeo virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Phasi Charoen-like phasivirus', 'Ochlerotatus scapularis flavivirus'] ,
('Australia', 23) : ['Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Phasi Charoen-like phasivirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Gabon', 2) : ['Mercadeo virus', 'Australian Anopheles totivirus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Tongilchon ohlsrhavirus', 'Zika virus', 'Menghai flavivirus'] ,
('Gabon', 8) : ['Pestivirus A', 'Alfalfa mosaic virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Merida virus', 'Lassa mammarenavirus', 'Aedes aegypti anphevirus', 'Mercadeo virus', 'Ohlsdorf ohlsrhavirus', 'Gata virus', 'Culex ohlsrhavirus', 'West Nile virus', 'Chinese rice-field eel rhabdovirus', 'Grenada mosquito rhabdovirus 1', 'Zika virus', 'Yata ephemerovirus', 'Perinet vesiculovirus', 'New Kent County virus'] ,
('Philippines', 24) : ['Wuhan Mosquito Virus 6', 'Liao ning virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Prunus necrotic ringspot virus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'West Nile virus', 'Mercadeo virus', 'Quang Binh virus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('Mexico', 19) : ['Liao ning virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Menghai flavivirus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Zika virus', 'Xishuangbanna aedes flavivirus'] ,
('South-Africa', 11) : ['Lassa mammarenavirus', 'Prunus necrotic ringspot virus', 'Alfalfa mosaic virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Radi vesiculovirus', 'Tongilchon ohlsrhavirus', 'Mercadeo virus', 'Zika virus', 'Phasi Charoen-like phasivirus'] ,
('Mexico', 21) : ['Lassa mammarenavirus', 'Liao ning virus', 'Merida virus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Menghai flavivirus', 'Avian orthoavulavirus 1', 'Mercadeo virus', 'Zika virus', 'Kamiti River virus', 'Pestivirus A'] ,
('Thailand', 18) : ['Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'West Nile virus', 'Pestivirus A', 'Dengue virus', 'Lassa mammarenavirus', 'Zika virus', 'Xingshan alphanemrhavirus'] ,
('French-Polynesia', 3) : ['Lassa mammarenavirus', 'Ochlerotatus caspius flavivirus', 'Ohlsdorf ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex flavivirus', 'Cell fusing agent virus', 'Zika virus', 'Riverside virus 1'] ,
('Philippines', 11) : ['Menghai flavivirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Mercadeo virus', 'Lassa mammarenavirus', 'Zika virus', 'Avian orthoavulavirus 1', 'Ochlerotatus scapularis flavivirus', 'West Nile virus', 'Ohlsdorf ohlsrhavirus'] ,
('Brazil', 8) : ['Liao ning virus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Merida virus', 'Aedes anphevirus', 'Alfalfa mosaic virus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Serbia mononega-like virus 1'] ,
('Brazil', 3) : ['Lassa mammarenavirus', 'Ohlsdorf ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Aedes anphevirus', 'Alfalfa mosaic virus', 'Culex tritaeniorhynchus rhabdovirus', 'Merida virus', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Phasi Charoen-like phasivirus', 'Fitzroy Crossing toti-like virus 1', 'Zika virus'] ,
('Mexico', 14) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Phasi Charoen-like phasivirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Liao ning virus', 'Aedes anphevirus'] ,
('South-Africa', 13) : ['Lassa mammarenavirus', 'Aedes aegypti anphevirus', 'Alfalfa mosaic virus', 'Aedes anphevirus', 'Tongilchon ohlsrhavirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Mercadeo virus', 'Zika virus', 'Pestivirus A'] ,
('Thailand', 10) : ['Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Lassa mammarenavirus', 'Mercadeo virus', 'Cell fusing agent virus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Menghai flavivirus'] ,
('Brazil', 11) : ['Lassa mammarenavirus', 'Pestivirus A', 'Tongilchon ohlsrhavirus', 'Liao ning virus', 'Aedes anphevirus', 'Zika virus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Kamiti River virus', 'Australian Anopheles totivirus'] ,
('Mexico', 17) : ['Lassa mammarenavirus', 'Liao ning virus', 'Aedes anphevirus', 'West Nile virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Mercadeo virus', 'Australian Anopheles totivirus'] ,
('USA', 3) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Mercadeo virus', 'Pestivirus A', 'Zika virus', 'Cell fusing agent virus', 'Avian orthoavulavirus 1'] ,
('Australia', 5) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Zika virus'] ,
('USA', 8) : ['Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Liao ning virus', 'West Nile virus', 'Ochlerotatus scapularis flavivirus', 'Lassa mammarenavirus', 'Zika virus', 'Avian orthoavulavirus 1', 'Pestivirus A'] ,
('USA', 14) : ['Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Pestivirus A', 'Ohlsdorf ohlsrhavirus', 'Riverside virus 1', 'Menghai flavivirus', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Sabethes flavivirus', 'Zika virus', 'Phasi Charoen-like phasivirus', 'Avian orthoavulavirus 1'] ,
('Thailand', 16) : ['Lassa mammarenavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Radi vesiculovirus', 'Ochlerotatus scapularis flavivirus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Pestivirus A'] ,
('Thailand', 8) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Wuhan Mosquito Virus 6', 'Menghai flavivirus', 'Aedes anphevirus', 'Liao ning virus', 'Pestivirus A', 'Australian Anopheles totivirus', 'Zika virus'] ,
('Angola', 18) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Zika virus', 'Ohlsdorf ohlsrhavirus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Aedes aegypti anphevirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Fitzroy Crossing toti-like virus 1', 'Cell fusing agent virus', 'Dengue virus'] ,
('Brazil', 10) : ['Lassa mammarenavirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Wuhan Mosquito Virus 6', 'Merida virus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Zika virus', 'Australian Anopheles totivirus', 'Ochlerotatus scapularis flavivirus'] ,
('Brazil', 14) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Kamiti River virus', 'Zika virus'] ,
('Angola', 19) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Aedes anphevirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Ochlerotatus caspius flavivirus', 'Culex ohlsrhavirus', 'Zika virus', 'Culiseta flavivirus', 'Xishuangbanna aedes flavivirus', 'Ochlerotatus scapularis flavivirus'] ,
('Thailand', 20) : ['Lassa mammarenavirus', 'Culiseta flavivirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Mexico', 23) : ['Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Phasi Charoen-like phasivirus', 'Ochlerotatus scapularis flavivirus', 'Menghai flavivirus', 'Lassa mammarenavirus', 'Culex flavivirus', 'Alfalfa mosaic virus', 'Xishuangbanna aedes flavivirus', 'Pestivirus A'] ,
('South-Africa', 10) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Ochlerotatus caspius flavivirus', 'Prunus necrotic ringspot virus', 'Wuhan Mosquito Virus 6', 'Radi vesiculovirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Phasi Charoen-like phasivirus', 'Culex ohlsrhavirus', 'Riverside virus 1', 'Zika virus', 'Culex rhabdovirus'] ,
('Angola', 2) : ['Ochlerotatus scapularis flavivirus', 'Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes aegypti anphevirus', 'Zika virus', 'Phasi Charoen-like phasivirus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus'] ,
('Vietnam', 11) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Liao ning virus', 'Tongilchon ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Kamiti River virus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Pestivirus A', 'West Nile virus'] ,
('Australia', 8) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Alfalfa mosaic virus', 'Aedes anphevirus', 'Liao ning virus', 'Cell fusing agent virus', 'Zika virus', 'Pestivirus A'] ,
('South-Africa', 5) : ['Lassa mammarenavirus', 'Tongilchon ohlsrhavirus', 'Alfalfa mosaic virus', 'Ohlsdorf ohlsrhavirus', 'Culex mosquito virus 2', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Riverside virus 1', 'Zika virus', 'Parramatta River virus', 'Avian orthoavulavirus 1'] ,
('French-Polynesia', 4) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes aegypti anphevirus', 'Liao ning virus', 'Tongilchon ohlsrhavirus', 'Zika virus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Japanese encephalitis virus', 'Cell fusing agent virus', 'Yongjia ledantevirus', 'Xishuangbanna aedes flavivirus'] ,
('Vietnam', 6) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Kamiti River virus', 'Zika virus', 'Ochlerotatus scapularis flavivirus', 'Pestivirus A'] ,
('Vietnam', 20) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Liao ning virus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Kamiti River virus', 'Zika virus'] ,
('Thailand', 2) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Pestivirus A'] ,
('Angola', 11) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes aegypti anphevirus', 'Tongilchon ohlsrhavirus', 'Zika virus', 'Ohlsdorf ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus'] ,
('Brazil', 12) : ['Lassa mammarenavirus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Alfalfa mosaic virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Zika virus', 'Fitzroy Crossing toti-like virus 1'] ,
('Philippines', 9) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'West Nile virus', 'Zika virus', 'Xincheng anphevirus'] ,
('South-Africa', 2) : ['Alfalfa mosaic virus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Culex mosquito virus 2', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Zika virus', 'Aedes anphevirus', 'Phasi Charoen-like phasivirus'] ,
('Philippines', 3) : ['Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Mercadeo virus', 'Lassa mammarenavirus', 'Zika virus'] ,
('Gabon', 23) : ['Lassa mammarenavirus', 'Culiseta flavivirus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus'] ,
('French-Polynesia', 7) : ['Tongilchon ohlsrhavirus', 'Lassa mammarenavirus', 'Cell fusing agent virus', 'Liao ning virus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Mercadeo virus', 'Xishuangbanna aedes flavivirus', 'Alfalfa mosaic virus', 'Zika virus', 'Phasi Charoen-like phasivirus', 'Island sawgrhavirus'] ,
('Mexico', 7) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Phasi Charoen-like phasivirus', 'La Tina virus'] ,
('Angola', 3) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Ohlsdorf ohlsrhavirus', 'Alfalfa mosaic virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Zika virus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Riverside virus 1'] ,
('French-Polynesia', 13) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Tongilchon ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Xishuangbanna aedes flavivirus', 'Ochlerotatus scapularis flavivirus', 'Avian orthoavulavirus 1', 'Alfalfa mosaic virus'] ,
('Vietnam', 7) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus'] ,
('Mexico', 4) : ['Culex ohlsrhavirus', 'Mercadeo virus', 'Aedes anphevirus', 'Lassa mammarenavirus', 'Liao ning virus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'Pestivirus A', 'Avian orthoavulavirus 1'] ,
('Philippines', 13) : ['Merida virus', 'Mercadeo virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Pestivirus A', 'West Nile virus', 'Zika virus', 'Menghai flavivirus'] ,
('South-Africa', 4) : ['Alfalfa mosaic virus', 'Prunus necrotic ringspot virus', 'Lassa mammarenavirus', 'Ohlsdorf ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Zika virus', 'Avian orthoavulavirus 1', 'Cell fusing agent virus'] ,
('Vietnam', 10) : ['Cell fusing agent virus', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Nienokoue virus', 'Prunus necrotic ringspot virus', 'Zika virus', 'Kamiti River virus'] ,
('Australia', 19) : ['Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Aedes anphevirus', 'Xishuangbanna aedes flavivirus', 'Lassa mammarenavirus', 'Zika virus'] ,
('Thailand', 19) : ['Mercadeo virus', 'Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Cell fusing agent virus', 'Xishuangbanna aedes flavivirus', 'Australian Anopheles totivirus', 'Zika virus'] ,
('Gabon', 1) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Merida virus', 'Culex tritaeniorhynchus rhabdovirus', 'Wuhan Mosquito Virus 6', 'Ochlerotatus scapularis flavivirus', 'Ohlsdorf ohlsrhavirus', 'Zika virus'] ,
('South-Africa', 7) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Culex mosquito virus 2', 'Radi vesiculovirus', 'Wuhan Mosquito Virus 6', 'Pestivirus A', 'Zika virus', 'Prunus necrotic ringspot virus'] ,
('Gabon', 15) : ['Pestivirus A', 'Tongilchon ohlsrhavirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Lassa mammarenavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Phasi Charoen-like phasivirus', 'Ohlsdorf ohlsrhavirus', 'Alfalfa mosaic virus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Liao ning virus', 'Culiseta flavivirus', 'Culex mosquito virus 2', 'Culex ohlsrhavirus', 'Mac Peak virus'] ,
('South-Africa', 3) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Culex tritaeniorhynchus rhabdovirus', 'Prunus necrotic ringspot virus', 'Zika virus', 'Alfalfa mosaic virus', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Cell fusing agent virus', 'Phasi Charoen-like phasivirus'] ,
('French-Polynesia', 19) : ['Cell fusing agent virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Tongilchon ohlsrhavirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Culex flavivirus', 'Australian Anopheles totivirus', 'Fitzroy Crossing toti-like virus 1', 'Ochlerotatus scapularis flavivirus', 'Barkedji virus', 'Zika virus'] ,
('Australia', 2) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Tongilchon ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Cell fusing agent virus', 'Zika virus'] ,
('Vietnam', 12) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Mercadeo virus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Ohlsdorf ohlsrhavirus', 'Zika virus', 'Pestivirus A', 'Australian Anopheles totivirus', 'Xishuangbanna aedes flavivirus'] ,
('Mexico', 15) : ['Mercadeo virus', 'Liao ning virus', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Phasi Charoen-like phasivirus', 'Australian Anopheles totivirus'] ,
('South-Africa', 15) : ['Lassa mammarenavirus', 'Prunus necrotic ringspot virus', 'Australian Anopheles totivirus', 'Tongilchon ohlsrhavirus', 'Alfalfa mosaic virus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Fitzroy Crossing toti-like virus 1', 'Aedes anphevirus', 'Phasi Charoen-like phasivirus', 'Zika virus', 'Ochlerotatus scapularis flavivirus'] ,
('Vietnam', 13) : ['Cell fusing agent virus', 'Culex ohlsrhavirus', 'Lassa mammarenavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Avian orthoavulavirus 1', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Menghai flavivirus'] ,
('South-Africa', 8) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Prunus necrotic ringspot virus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes aegypti anphevirus', 'Culiseta flavivirus', 'Culex ohlsrhavirus', 'Zika virus', 'West Nile virus', 'Culex flavivirus'] ,
('USA', 13) : ['Liao ning virus', 'Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Culex ohlsrhavirus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus'] ,
('Thailand', 4) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Culex ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Australian Anopheles totivirus', 'Pestivirus A', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Australia', 24) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Liao ning virus', 'Culiseta flavivirus', 'Ohlsdorf ohlsrhavirus', 'Alfalfa mosaic virus', 'Australian Anopheles totivirus', 'Pestivirus A', 'Zika virus'] ,
('Philippines', 19) : ['Liao ning virus', 'Wuhan Mosquito Virus 6', 'Lassa mammarenavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'West Nile virus', 'Pestivirus A', 'Riverside virus 1', 'Mercadeo virus', 'Zika virus', 'Menghai flavivirus', 'Xingshan alphanemrhavirus'] ,
('Mexico', 13) : ['Lassa mammarenavirus', 'Ochlerotatus scapularis flavivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Avian orthoavulavirus 1'] ,
('Gabon', 4) : ['Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Culex ohlsrhavirus', 'Culiseta flavivirus', 'Culex tritaeniorhynchus rhabdovirus'] ,
('Australia', 22) : ['Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Liao ning virus', 'Tongilchon ohlsrhavirus', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Mercadeo virus', 'Xishuangbanna aedes flavivirus', 'Riverside virus 1', 'Lassa mammarenavirus', 'Zika virus', 'Aedes albopictus anphevirus', 'Pestivirus A'] ,
('Australia', 21) : ['Ochlerotatus scapularis flavivirus', 'Lassa mammarenavirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Australian Anopheles totivirus', 'Zika virus'] ,
('Australia', 13) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Liao ning virus', 'Alfalfa mosaic virus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Avian orthoavulavirus 1', 'Zika virus', 'La Tina virus'] ,
('Thailand', 12) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'West Nile virus', 'Australian Anopheles totivirus', 'Kamiti River virus', 'Zika virus', 'Xishuangbanna aedes flavivirus'] ,
('Vietnam', 14) : ['Lassa mammarenavirus', 'Liao ning virus', 'Parramatta River virus', 'Cell fusing agent virus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Phasi Charoen-like phasivirus', 'Zika virus'] ,
('Philippines', 7) : ['Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Liao ning virus', 'Pestivirus A', 'Cell fusing agent virus', 'Avian orthoavulavirus 1', 'Zika virus'] ,
('French-Polynesia', 16) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Aedes aegypti anphevirus', 'Ohlsdorf ohlsrhavirus', 'Riverside virus 1', 'Culex ohlsrhavirus', 'Liao ning virus', 'West Nile virus', 'Alfalfa mosaic virus', 'Zika virus', 'Xishuangbanna aedes flavivirus', 'Pestivirus A'] ,
('USA', 6) : ['Merida virus', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Mercadeo virus', 'Liao ning virus', 'Avian orthoavulavirus 1', 'Pestivirus A', 'Culex ohlsrhavirus'] ,
('Philippines', 18) : ['Liao ning virus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Cell fusing agent virus', 'Tongilchon ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Pestivirus A', 'Mercadeo virus', 'West Nile virus', 'Zika virus', 'Ochlerotatus scapularis flavivirus'] ,
('Thailand', 6) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Tongilchon ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Pestivirus A', 'Zika virus'] ,
('South-Africa', 6) : ['Alfalfa mosaic virus', 'Lassa mammarenavirus', 'Prunus necrotic ringspot virus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Culex mosquito virus 2', 'Culex tritaeniorhynchus rhabdovirus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'Aedes aegypti anphevirus', 'Zika virus', 'Aedes anphevirus', 'Riverside virus 1', 'Yinshui bat virus', 'Pestivirus A'] ,
('Thailand', 3) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Parramatta River virus', 'Australian Anopheles totivirus', 'Zika virus'] ,
('Argentina', 31) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Ohlsdorf ohlsrhavirus', 'Zika virus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Riverside virus 1', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Fitzroy Crossing toti-like virus 1', 'Ochlerotatus scapularis flavivirus', 'Xishuangbanna aedes flavivirus', 'Inhangapi virus', 'Prunus necrotic ringspot virus'] ,
('Thailand', 9) : ['Aedes anphevirus', 'Cell fusing agent virus', 'Liao ning virus', 'Menghai flavivirus', 'Mercadeo virus', 'Lassa mammarenavirus', 'Xishuangbanna aedes flavivirus', 'Zika virus'] ,
('Brazil', 18) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Ochlerotatus scapularis flavivirus', 'Xishuangbanna aedes flavivirus', 'Avian orthoavulavirus 1', 'Australian Anopheles totivirus', 'Liao ning virus', 'Serbia mononega-like virus 1'] ,
('USA', 21) : ['Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Aedes anphevirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Pestivirus A', 'Culex mosquito virus 2', 'Zika virus', 'Parramatta River virus'] ,
('French-Polynesia', 17) : ['Lassa mammarenavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Australian Anopheles totivirus', 'Tongilchon ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Zika virus', 'Pestivirus A', 'Alfalfa mosaic virus'] ,
('Angola', 6) : ['Lassa mammarenavirus', 'Alfalfa mosaic virus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Riverside virus 1', 'Zika virus', 'Cell fusing agent virus', 'Ohlsdorf ohlsrhavirus', 'Radi vesiculovirus', 'Liao ning virus', 'Mercadeo virus', 'Pestivirus A', 'Avian orthoavulavirus 1'] ,
('Gabon', 13) : ['Lassa mammarenavirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Phasi Charoen-like phasivirus', 'Aedes anphevirus', 'Ohlsdorf ohlsrhavirus', 'Cell fusing agent virus', 'West Nile virus', 'Yinshui bat virus', 'Culex rhabdovirus', 'Guadeloupe Culex rhabdovirus', 'Kamiti River virus', 'Avian orthoavulavirus 1'] ,
('French-Polynesia', 8) : ['Lassa mammarenavirus', 'Cell fusing agent virus', 'Wuhan Mosquito Virus 6', 'Ohlsdorf ohlsrhavirus', 'Radi vesiculovirus', 'Ochlerotatus scapularis flavivirus', 'Australian Anopheles totivirus', 'Aedes aegypti anphevirus', 'Liao ning virus', 'Xishuangbanna aedes flavivirus', 'Alfalfa mosaic virus', 'Zika virus', 'Prunus necrotic ringspot virus'] ,
('Vietnam', 5) : ['Lassa mammarenavirus', 'Mercadeo virus', 'Cell fusing agent virus', 'Australian Anopheles totivirus', 'Aedes anphevirus', 'Culex tritaeniorhynchus rhabdovirus', 'Liao ning virus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Zika virus'] ,
('Angola', 20) : ['Ohlsdorf ohlsrhavirus', 'Liao ning virus', 'Alfalfa mosaic virus', 'Wuhan Mosquito Virus 6', 'Lassa mammarenavirus', 'Aedes anphevirus', 'Culex ohlsrhavirus', 'Zika virus', 'Radi vesiculovirus', 'Fitzroy Crossing toti-like virus 1', 'Dengue virus', 'Australian Anopheles totivirus', 'Avian orthoavulavirus 1'] ,
('Angola', 4) : ['Lassa mammarenavirus', 'Liao ning virus', 'Culex Flavi-like virus', 'Ohlsdorf ohlsrhavirus', 'Aedes anphevirus', 'Aedes aegypti anphevirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Culex tritaeniorhynchus rhabdovirus', 'Wuhan Mosquito Virus 6', 'Australian Anopheles totivirus', 'Zika virus', 'Pestivirus A'] ,
('Gabon', 22) : ['Lassa mammarenavirus', 'Menghai flavivirus', 'Australian Anopheles totivirus', 'Alfalfa mosaic virus', 'Tongilchon ohlsrhavirus', 'West Nile virus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Phasi Charoen-like phasivirus', 'Wuhan Mosquito Virus 6', 'Liao ning virus', 'Mercadeo virus', 'Pestivirus A', 'Zika virus', 'Aedes anphevirus'] ,
('Australia', 15) : ['Lassa mammarenavirus', 'Wuhan Mosquito Virus 6', 'Aedes anphevirus', 'Mercadeo virus', 'Australian Anopheles totivirus', 'Culex ohlsrhavirus', 'Cell fusing agent virus', 'Liao ning virus', 'Riverside virus 1', 'Xishuangbanna aedes flavivirus', 'Zika virus', 'Avian orthoavulavirus 1'] ,
('French-Polynesia', 1) : ['Lassa mammarenavirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Culex tritaeniorhynchus rhabdovirus', 'Tongilchon ohlsrhavirus', 'Liao ning virus', 'Cell fusing agent virus', 'Culex ohlsrhavirus', 'Wuhan Mosquito Virus 6', 'Xishuangbanna aedes flavivirus', 'Ochlerotatus scapularis flavivirus', 'Xingshan alphanemrhavirus', 'Alfalfa mosaic virus', 'Zika virus'] ,
('Mexico', 6) : ['Liao ning virus', 'Wuhan Mosquito Virus 6', 'Lassa mammarenavirus', 'Culex ohlsrhavirus', 'Avian orthoavulavirus 1', 'Aedes anphevirus', 'Menghai flavivirus', 'Australian Anopheles totivirus', 'Cell fusing agent virus', 'Zika virus', 'Pestivirus A', 'Mercadeo virus', 'Xishuangbanna aedes flavivirus'] ,
('Gabon', 20) : ['Lassa mammarenavirus', 'Culiseta flavivirus', 'Aedes anphevirus', 'Australian Anopheles totivirus', 'Wuhan Mosquito Virus 6', 'Tongilchon ohlsrhavirus', 'Pestivirus A', 'Ohlsdorf ohlsrhavirus', 'Xishuangbanna aedes flavivirus', 'Xingshan alphanemrhavirus', 'Zika virus', 'Ochlerotatus scapularis flavivirus']
}

all_viruses = ['Aedes aegypti anphevirus', 'Aedes anphevirus', 'American plum line pattern virus', 'Australian Anopheles totivirus',
 'Cell fusing agent virus', 'Culex Flavi-like virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culiseta flavivirus',
  'Dengue virus', 'Kamiti River virus', 'Lassa mammarenavirus', 'Liao ning virus', 'Mercadeo virus', 'Merida virus',
   'Ochlerotatus caspius flavivirus', 'Ochlerotatus scapularis flavivirus', 'Ohlsdorf ohlsrhavirus', 'Pestivirus A', 'Radi vesiculovirus',
   'Tongilchon ohlsrhavirus', 'West Nile virus', 'Wuhan Mosquito Virus 6', 'Xishuangbanna aedes flavivirus']

vFam = {'Arenaviridae': ['Lassa mammarenavirus'], 'Bromoviridae': ['Alfalfa mosaic virus', 'American plum line pattern virus', 'Prunus necrotic ringspot virus'], 'Flaviviridae': ['Aedes flavivirus', 'Barkedji virus', 'Cell fusing agent virus', 'Culex Flavi-like virus', 'Culex flavivirus', 'Culiseta flavivirus', 'Dengue virus', 'Japanese encephalitis virus', 'Kamiti River virus', 'La Tina virus', 'Mac Peak virus', 'Mansonia flavivirus', 'Menghai flavivirus', 'Mercadeo virus', 'Nienokoue virus', 'Ochlerotatus caspius flavivirus', 'Ochlerotatus scapularis flavivirus', 'Parramatta River virus', 'Pestivirus A', 'Quang Binh virus', 'Sabethes flavivirus', 'Tick-borne encephalitis virus', 'West Nile virus', 'Xishuangbanna aedes flavivirus', 'Zika virus'], 'Mononegavirales': ['Serbia mononega-like virus 1'], 'Nodaviridae': ['Culex mosquito virus 2'], 'Orthomyxoviridae': ['Wuhan Mosquito Virus 6'], 'Paramyxoviridae': ['Avian orthoavulavirus 1'], 'Phenuiviridae': ['Phasi Charoen-like phasivirus'], 'Pneumoviridae': ['Human orthopneumovirus'], 'Rhabdoviridae': ['Chinese rice-field eel rhabdovirus', 'Culex ohlsrhavirus', 'Culex rhabdovirus', 'Culex tritaeniorhynchus rhabdovirus', 'Drosophila melanogaster sigmavirus', 'Gata virus', 'Grenada mosquito rhabdovirus 1', 'Guadeloupe Culex rhabdovirus', 'Inhangapi virus', 'Island sawgrhavirus', 'Malpais Spring vesiculovirus', 'Merida virus', 'New Kent County virus', 'Ohlsdorf ohlsrhavirus', 'Perinet vesiculovirus', 'Radi vesiculovirus', 'Riverside virus 1', 'Tongilchon ohlsrhavirus', 'Xingshan alphanemrhavirus', 'Yata ephemerovirus', 'Yinshui bat virus', 'Yongjia ledantevirus'], 'Sedoreovirinae': ['Liao ning virus'], 'Totiviridae': ['Australian Anopheles totivirus', 'Fitzroy Crossing toti-like virus 1'], 'Xinmoviridae': ['Aedes aegypti anphevirus', 'Aedes albopictus anphevirus', 'Aedes anphevirus', 'Xincheng anphevirus']}


all_regions = {
    'Angola': [1,2,3,4,5,6,8,9,10,11,12,14,16,17,18,19,20],
    'Argentina': [2,3,8,11,31],
    'Australia': [1,2,5,7,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24],
    'Brazil': [1,2,3,4,6,7,8,10,11,12,13,14,15,18,19,23,24],
    'French-Polynesia': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,10],
    'Gabon': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24],
    'Mexico': [1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
    'Philippines': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24],
    'South-Africa': [2,3,4,5,6,7,8,9,10,11,12,13,14,15],
    'Thailand': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
    'USA': [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,18,19,20,21,22,23,24],
    'Vietnam': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
}

label2Specimen = {('Mexico', 1): 'Amacuzac-Mexico-1.LIN210A1618', ('Mexico', 10): 'Amacuzac-Mexico-10.LIN210A1625', ('Mexico', 11): 'Amacuzac-Mexico-11.LIN210A1626', ('Mexico', 13): 'Amacuzac-Mexico-13.LIN210A1627', ('Mexico', 14): 'Amacuzac-Mexico-14.LIN210A1628', ('Mexico', 16): 'Amacuzac-Mexico-16.LIN210A1629', ('Mexico', 17): 'Amacuzac-Mexico-17.LIN210A1630', ('Mexico', 18): 'Amacuzac-Mexico-18.LIN210A1631', ('Mexico', 19): 'Amacuzac-Mexico-19.LIN210A1632', ('Mexico', 20): 'Amacuzac-Mexico-20.LIN210A1633', ('Mexico', 22): 'Amacuzac-Mexico-22.LIN210A1634', ('Mexico', 23): 'Amacuzac-Mexico-23.LIN210A1635', ('Mexico', 24): 'Amacuzac-Mexico-24.LIN210A1636', ('Mexico', 3): 'Amacuzac-Mexico-3.LIN210A1619', ('Mexico', 4): 'Amacuzac-Mexico-4.LIN210A1620', ('Mexico', 5): 'Amacuzac-Mexico-5.LIN210A1621', ('Mexico', 6): 'Amacuzac-Mexico-6.LIN210A1622', ('Mexico', 7): 'Amacuzac-Mexico-7.LIN210A1623', ('Mexico', 8): 'Amacuzac-Mexico-8.LIN210A1624', ('Mexico', 15): 'Amacuzac_Mexico-_15.2-116629', ('Mexico', 2): 'Amacuzac_Mexico-_2.2-116616', ('Mexico', 21): 'Amacuzac_Mexico-_21.2-116635', ('Thailand', 1): 'Bangkok_Thailand_01.LIN210A1679', ('Thailand', 2): 'Bangkok_Thailand_02.LIN210A1680', ('Thailand', 3): 'Bangkok_Thailand_03.LIN210A1681', ('Thailand', 4): 'Bangkok_Thailand_04.LIN210A1682', ('Thailand', 5): 'Bangkok_Thailand_05.LIN210A1683', ('Thailand', 6): 'Bangkok_Thailand_06.LIN210A1684', ('Thailand', 7): 'Bangkok_Thailand_07.LIN210A1685', ('Thailand', 8): 'Bangkok_Thailand_08.LIN210A1686', ('Thailand', 9): 'Bangkok_Thailand_09.LIN210A1687', ('Thailand', 10): 'Bangkok_Thailand_10.LIN210A1688', ('Thailand', 11): 'Bangkok_Thailand_11.LIN210A1689', ('Thailand', 12): 'Bangkok_Thailand_12.LIN210A1690', ('Thailand', 13): 'Bangkok_Thailand_13.LIN210A1691', ('Thailand', 14): 'Bangkok_Thailand_14.LIN210A1692', ('Thailand', 15): 'Bangkok_Thailand_15.LIN210A1693', ('Thailand', 16): 'Bangkok_Thailand_16.LIN210A1694', ('Thailand', 17): 'Bangkok_Thailand_17.LIN210A1695', ('Thailand', 18): 'Bangkok_Thailand_18.LIN210A1696', ('Thailand', 19): 'Bangkok_Thailand_19.LIN210A1697', ('Thailand', 20): 'Bangkok_Thailand_20.LIN210A1698', ('Australia', 1): 'Cairns-Australia-1.LIN210A2151', ('Australia', 10): 'Cairns-Australia-10.LIN210A2160', ('Australia', 11): 'Cairns-Australia-11.LIN210A2161', ('Australia', 12): 'Cairns-Australia-12.LIN210A2162', ('Australia', 13): 'Cairns-Australia-13.LIN210A2163', ('Australia', 15): 'Cairns-Australia-15.LIN210A2165', ('Australia', 16): 'Cairns-Australia-16.LIN210A2166', ('Australia', 17): 'Cairns-Australia-17.LIN210A2167', ('Australia', 18): 'Cairns-Australia-18.LIN210A2168', ('Australia', 19): 'Cairns-Australia-19.LIN210A2169', ('Australia', 2): 'Cairns-Australia-2.LIN210A2152', ('Australia', 20): 'Cairns-Australia-20.LIN210A2170', ('Australia', 21): 'Cairns-Australia-21.LIN210A2171', ('Australia', 22): 'Cairns-Australia-22.LIN210A2172', ('Australia', 23): 'Cairns-Australia-23.LIN210A2173', ('Australia', 24): 'Cairns-Australia-24.LIN210A2174', ('Australia', 5): 'Cairns-Australia-5.LIN210A2155', ('Australia', 7): 'Cairns-Australia-7.LIN210A2157', ('Australia', 8): 'Cairns-Australia-8.LIN210A2158', ('Philippines', 1): 'Cebu_City_Philippines-_1.2-116713', ('Philippines', 10): 'Cebu_City_Philippines-_10.2-116722', ('Philippines', 11): 'Cebu_City_Philippines-_11.2-116723', ('Philippines', 12): 'Cebu_City_Philippines-_12.2-116724', ('Philippines', 13): 'Cebu_City_Philippines-_13.2-116725', ('Philippines', 14): 'Cebu_City_Philippines-_14.2-116726', ('Philippines', 15): 'Cebu_City_Philippines-_15.2-116727', ('Philippines', 16): 'Cebu_City_Philippines-_16.2-116728', ('Philippines', 18): 'Cebu_City_Philippines-_18.2-116730', ('Philippines', 19): 'Cebu_City_Philippines-_19.2-116731', ('Philippines', 2): 'Cebu_City_Philippines-_2.2-116714', ('Philippines', 20): 'Cebu_City_Philippines-_20.2-116732', ('Philippines', 21): 'Cebu_City_Philippines-_21.2-116733', ('Philippines', 22): 'Cebu_City_Philippines-_22.2-116734', ('Philippines', 23): 'Cebu_City_Philippines-_23.2-116735', ('Philippines', 24): 'Cebu_City_Philippines-_24.2-116736', ('Philippines', 3): 'Cebu_City_Philippines-_3.2-116715', ('Philippines', 4): 'Cebu_City_Philippines-_4.2-116716', ('Philippines', 5): 'Cebu_City_Philippines-_5.2-116717', ('Philippines', 6): 'Cebu_City_Philippines-_6.2-116718', ('Philippines', 7): 'Cebu_City_Philippines-_7.2-116719', ('Philippines', 8): 'Cebu_City_Philippines-_8.2-116720', ('Philippines', 9): 'Cebu_City_Philippines-_9.2-116721', ('Angola', 1): 'Cuanda_Angola_01.LIN210A1719', ('Angola', 2): 'Cuanda_Angola_02.LIN210A1720', ('Angola', 3): 'Cuanda_Angola_03.LIN210A1721', ('Angola', 4): 'Cuanda_Angola_04.LIN210A1722', ('Angola', 5): 'Cuanda_Angola_05.LIN210A1723', ('Angola', 6): 'Cuanda_Angola_06.LIN210A1724', ('Angola', 8): 'Cuanda_Angola_08.LIN210A1726', ('Angola', 9): 'Cuanda_Angola_09.LIN210A1727', ('Angola', 10): 'Cuanda_Angola_10.LIN210A1728', ('Angola', 11): 'Cuanda_Angola_11.LIN210A1729', ('Angola', 12): 'Cuanda_Angola_12.LIN210A1730', ('Angola', 14): 'Cuanda_Angola_14.LIN210A1732', ('Angola', 16): 'Cuanda_Angola_16.LIN210A1734', ('Angola', 17): 'Cuanda_Angola_17.LIN210A1735', ('Angola', 18): 'Cuanda_Angola_18.LIN210A1736', ('Angola', 19): 'Cuanda_Angola_19.LIN210A1737', ('Angola', 20): 'Cuanda_Angola_20.LIN210A1738', ('Argentina', 11): 'El_Dorado_Argentina_U_11.LIN210A2770', ('Argentina', 2): 'El_Dorado_Argentina_U_2.LIN210A2698', ('Argentina', 3): 'El_Dorado_Argentina_U_3.LIN210A2706', ('Argentina', 8): 'El_Dorado_Argentina_U_8.LIN210A2746', ('Argentina', 31): 'El_Dorado_US_U_31', ('Vietnam', 1): 'HoChiMin_Vietnam_01.LIN210A1659', ('Vietnam', 2): 'HoChiMin_Vietnam_02.LIN210A1660', ('Vietnam', 3): 'HoChiMin_Vietnam_03.LIN210A1661', ('Vietnam', 4): 'HoChiMin_Vietnam_04.LIN210A1662', ('Vietnam', 5): 'HoChiMin_Vietnam_05.LIN210A1663', ('Vietnam', 6): 'HoChiMin_Vietnam_06.LIN210A1664', ('Vietnam', 7): 'HoChiMin_Vietnam_07.LIN210A1665', ('Vietnam', 8): 'HoChiMin_Vietnam_08.LIN210A1666', ('Vietnam', 9): 'HoChiMin_Vietnam_09.LIN210A1667', ('Vietnam', 10): 'HoChiMin_Vietnam_10.LIN210A1668', ('Vietnam', 11): 'HoChiMin_Vietnam_11.LIN210A1669', ('Vietnam', 12): 'HoChiMin_Vietnam_12.LIN210A1670', ('Vietnam', 13): 'HoChiMin_Vietnam_13.LIN210A1671', ('Vietnam', 14): 'HoChiMin_Vietnam_14.LIN210A1672', ('Vietnam', 15): 'HoChiMin_Vietnam_15.LIN210A1673', ('Vietnam', 16): 'HoChiMin_Vietnam_16.LIN210A1674', ('Vietnam', 17): 'HoChiMin_Vietnam_17.LIN210A1675', ('Vietnam', 18): 'HoChiMin_Vietnam_18.LIN210A1676', ('Vietnam', 19): 'HoChiMin_Vietnam_19.LIN210A1677', ('Vietnam', 20): 'HoChiMin_Vietnam_20.LIN210A1678', ('Gabon', 1): 'La_Lope-Gabon-1.LIN210A1637', ('Gabon', 10): 'La_Lope-Gabon-10.LIN210A1646', ('Gabon', 11): 'La_Lope-Gabon-11.LIN210A1647', ('Gabon', 12): 'La_Lope-Gabon-12.LIN210A1648', ('Gabon', 13): 'La_Lope-Gabon-13.LIN210A1649', ('Gabon', 14): 'La_Lope-Gabon-14.LIN210A1650', ('Gabon', 15): 'La_Lope-Gabon-15.LIN210A1651', ('Gabon', 16): 'La_Lope-Gabon-16.LIN210A1652', ('Gabon', 17): 'La_Lope-Gabon-17.LIN210A1653', ('Gabon', 18): 'La_Lope-Gabon-18.LIN210A1654', ('Gabon', 2): 'La_Lope-Gabon-2.LIN210A1638', ('Gabon', 20): 'La_Lope-Gabon-20.LIN210A1655', ('Gabon', 23): 'La_Lope-Gabon-23.LIN210A1656', ('Gabon', 24): 'La_Lope-Gabon-24.LIN210A1658', ('Gabon', 3): 'La_Lope-Gabon-3.LIN210A1639', ('Gabon', 4): 'La_Lope-Gabon-4.LIN210A1640', ('Gabon', 5): 'La_Lope-Gabon-5.LIN210A1641', ('Gabon', 6): 'La_Lope-Gabon-6.LIN210A1642', ('Gabon', 7): 'La_Lope-Gabon-7.LIN210A1643', ('Gabon', 8): 'La_Lope-Gabon-8.LIN210A1644', ('Gabon', 9): 'La_Lope-Gabon-9.LIN210A1645', ('Gabon', 19): 'La_Lope_Gabon-_19.2-116707', ('Gabon', 22): 'La_Lope_Gabon-_22.2-116710', ('USA', 1): 'Maricopa_County_AZ-_1.2-116541', ('USA', 12): 'Maricopa_County_AZ-_12.2-116602', ('USA', 13): 'Maricopa_County_AZ-_13.2-116603', ('USA', 14): 'Maricopa_County_AZ-_14.2-116604', ('USA', 15): 'Maricopa_County_AZ-_15.2-116605', ('USA', 16): 'Maricopa_County_AZ-_16.2-116606', ('USA', 18): 'Maricopa_County_AZ-_18.2-116608', ('USA', 19): 'Maricopa_County_AZ-_19.2-116609', ('USA', 21): 'Maricopa_County_AZ-_21.2-116611', ('USA', 24): 'Maricopa_County_AZ-_24.2-116614', ('USA', 3): 'Maricopa_County_AZ-_3.2-116543', ('USA', 4): 'Maricopa_County_AZ-_4.2-116544', ('USA', 6): 'Maricopa_County_AZ-_6.2-116546', ('USA', 8): 'Maricopa_County_AZ-_8.2-116548', ('USA', 9): 'Maricopa_County_AZ-_9.2-116549', ('USA', 10): 'Maricopa_County_Arizona-USA-10.LIN210A1614', ('USA', 2): 'Maricopa_County_Arizona-USA-2.LIN210A1611', ('USA', 20): 'Maricopa_County_Arizona-USA-20.LIN210A1615', ('USA', 22): 'Maricopa_County_Arizona-USA-22.LIN210A1616', ('USA', 23): 'Maricopa_County_Arizona-USA-23.LIN210A1617', ('USA', 5): 'Maricopa_County_Arizona-USA-5.LIN210A1612', ('USA', 7): 'Maricopa_County_Arizona-USA-7.LIN210A1613', ('Brazil', 1): 'Paqueta-Brazil-1.LIN210A2127', ('Brazil', 10): 'Paqueta-Brazil-10.LIN210A2136', ('Brazil', 11): 'Paqueta-Brazil-11.LIN210A2137', ('Brazil', 12): 'Paqueta-Brazil-12.LIN210A2138', ('Brazil', 13): 'Paqueta-Brazil-13.LIN210A2139', ('Brazil', 14): 'Paqueta-Brazil-14.LIN210A2140', ('Brazil', 15): 'Paqueta-Brazil-15.LIN210A2141', ('Brazil', 18): 'Paqueta-Brazil-18.LIN210A2144', ('Brazil', 19): 'Paqueta-Brazil-19.LIN210A2145', ('Brazil', 2): 'Paqueta-Brazil-2.LIN210A2128', ('Brazil', 23): 'Paqueta-Brazil-23.LIN210A2149', ('Brazil', 24): 'Paqueta-Brazil-24.LIN210A2150', ('Brazil', 3): 'Paqueta-Brazil-3.LIN210A2129', ('Brazil', 4): 'Paqueta-Brazil-4.LIN210A2130', ('Brazil', 6): 'Paqueta-Brazil-6.LIN210A2132', ('Brazil', 7): 'Paqueta-Brazil-7.LIN210A2133', ('Brazil', 8): 'Paqueta-Brazil-8.LIN210A2134', ('South-Africa', 2): 'Skukusa_South_Africa_02.LIN210A1740', ('South-Africa', 3): 'Skukusa_South_Africa_03.LIN210A1741', ('South-Africa', 4): 'Skukusa_South_Africa_04.LIN210A1742', ('South-Africa', 5): 'Skukusa_South_Africa_05.LIN210A1743', ('South-Africa', 6): 'Skukusa_South_Africa_06.LIN210A1744', ('South-Africa', 7): 'Skukusa_South_Africa_07.LIN210A1745', ('South-Africa', 8): 'Skukusa_South_Africa_08.LIN210A1746', ('South-Africa', 9): 'Skukusa_South_Africa_09.LIN210A1747', ('South-Africa', 10): 'Skukusa_South_Africa_10.LIN210A1748', ('South-Africa', 11): 'Skukusa_South_Africa_11.LIN210A1749', ('South-Africa', 12): 'Skukusa_South_Africa_12.LIN210A1750', ('South-Africa', 13): 'Skukusa_South_Africa_13.LIN210A1751', ('South-Africa', 14): 'Skukusa_South_Africa_14.LIN210A1752', ('South-Africa', 15): 'Skukusa_South_Africa_15.LIN210A1753', ('French-Polynesia', 1): 'Tahiti_FrenchPolynesia_01.LIN210A1699', ('French-Polynesia', 2): 'Tahiti_FrenchPolynesia_02.LIN210A1700', ('French-Polynesia', 3): 'Tahiti_FrenchPolynesia_03.LIN210A1701', ('French-Polynesia', 4): 'Tahiti_FrenchPolynesia_04.LIN210A1702', ('French-Polynesia', 5): 'Tahiti_FrenchPolynesia_05.LIN210A1703', ('French-Polynesia', 6): 'Tahiti_FrenchPolynesia_06.LIN210A1704', ('French-Polynesia', 7): 'Tahiti_FrenchPolynesia_07.LIN210A1705', ('French-Polynesia', 8): 'Tahiti_FrenchPolynesia_08.LIN210A1706', ('French-Polynesia', 9): 'Tahiti_FrenchPolynesia_09.LIN210A1707', ('French-Polynesia', 10): 'Tahiti_FrenchPolynesia_10.LIN210A1708', ('French-Polynesia', 11): 'Tahiti_FrenchPolynesia_11.LIN210A1709', ('French-Polynesia', 12): 'Tahiti_FrenchPolynesia_12.LIN210A1710', ('French-Polynesia', 13): 'Tahiti_FrenchPolynesia_13.LIN210A1711', ('French-Polynesia', 14): 'Tahiti_FrenchPolynesia_14.LIN210A1712', ('French-Polynesia', 15): 'Tahiti_FrenchPolynesia_15.LIN210A1713', ('French-Polynesia', 16): 'Tahiti_FrenchPolynesia_16.LIN210A1714', ('French-Polynesia', 17): 'Tahiti_FrenchPolynesia_17.LIN210A1715', ('French-Polynesia', 18): 'Tahiti_FrenchPolynesia_18.LIN210A1716', ('French-Polynesia', 19): 'Tahiti_FrenchPolynesia_19.LIN210A1717', ('French-Polynesia', 20): 'Tahiti_FrenchPolynesia_20.LIN210A1718'}

regions = list(all_regions.keys())
regions.sort()

app = dash.Dash()

app.layout = html.Div(children=[
    html.H1(children='EVEs In Aedes aegypti Mosquitoes',style={'textAlign': 'center'}),

    html.Div(style={'textAlign': 'center'}, children=['Contig Diagrams for EVE in AAa genome',]),

    html.P("To enable more convenient appraisal of individual contigs, EVE produces diagrams of the viral hits in each contig (only the best hit). To see different contigs containing viral hits, select a specimen and virus.",
     style = {'textAlign': 'center'}),
    html.Br(),

    html.Div ([dcc.Dropdown( id='countries-dropdown',placeholder = 'Select a region', clearable = False, options=[{'label': k, 'value': k} for k in all_regions.keys()])], style = {'width':'20%', 'display': 'inline-block'},),

    html.Div([dcc.Dropdown(id='numbers-dropdown', placeholder = 'Select a num', clearable = False)], style = {'width':'20%', 'display': 'inline-block'},),

    html.Div([dcc.Dropdown(id='family-dropdown', options = [{'label': k, 'value': k} for k in vFam.keys()], clearable = False, placeholder = 'Select a family')], style = {'width': '30%', 'display': 'inline-block'}),

    html.Div([dcc.Dropdown(id='virus-dropdown', clearable = False, placeholder = 'Select a virus')], style = {'width':'30%', 'display': 'inline-block'}),

    html.Br(),html.Br(),

    dbc.Button('submit', id = 'submit', n_clicks = 0),

    html.P(id = 'result', style = {'display':'none'}),

    html.Br(), html.Br(),

    html.Div(id = 'contig-div', children = [html.P(id = 'choose-contig', style = {'font-size':14}),
            dcc.Input(id='contig-num', type ='number', placeholder = 'choose a contig', min=0, max = 5, value = 1),], style = {'display':'none'}),

    html.Div (id = 'match-div', children = [html.P(id = 'choose-matches', style = {'font-size':14}), dcc.Input(id='matches-num', type ='number', placeholder = 'Number of inv', min=0,  value = 0, )],
            style = {'display':'none'}),

    html.Div (id = 'inv-div', children =[html.P(id = 'choose-inv', style = {'font-size':14}), dcc.Input(id='inv-num', type ='number', placeholder = 'Number of inv', min=0, value = 0)],
             style = {'display':'none', "margin-left": "20px"}),

    html.Div (id = 'left-div', children =[html.P(id = 'choose-left', style = {'font-size':14}), dcc.Input(id='left-num', type ='number', placeholder = 'Number of inv', min=0, value = 0)],
             style = {'display':'none', "margin-left": "20px"}),

    html.Div (id = 'right-div', children =[html.P(id = 'choose-right', style = {'font-size':14}), dcc.Input(id='right-num', type ='number', placeholder = 'Number of inv', min=0, value = 0)],
             style = {'display':'none', "margin-left": "20px"}),

    html.Hr(),

    html.Div([html.Img(id = 'contigDiagram')], style = {'textAlign': 'center'}),

    html.Br(),

    html.Div(id = 'buttons', children = [
          dbc.Button('prev', id='left-scroll',n_clicks=0),
          dbc.Button('next', id = 'right-scroll', n_clicks = 0),],
          style = {'display':'none'}),

    html.P(id="download-seq", style = {'display': 'none'}),

    html.Button("Download Sequence", id="btn_txt", style = {'display': 'none'}),

    dcc.Download(id='download-text'),


])
###############################################################################

@app.callback(
    Output("download-text", "data"),
    Input("btn_txt", "n_clicks"),
    State('download-seq','children'),
    State('countries-dropdown', 'value'),
    State('numbers-dropdown', 'value'),
    prevent_initial_call=True,
)
def func(n_clicks, seq, country, num):
    contig = seq.split("|")[0]
    seqname = f"{country}_{num}_{contig}.fasta"
    return dict(content=seq, filename=seqname)

###############################################################################

@app.callback(
    Output('btn_txt', 'style'),
    Input('submit','n_clicks'),
    prevent_initial_call=True,
)
def hideDownload(n_clicks):
    return {'display' : 'inline-block'}

###############################################################################

@app.callback(
    Output('choose-matches', 'children'),
    Output('match-div', 'style'),
    Input('result', 'children'),
    State('matches-num', 'max'),
    prevent_initial_call = True
)
def set_matchNum_caption(submit, matches):
    return f"Choose number of matches (max: {matches})", {'display':'inline-block', "margin-left": "20px"}

###############################################################################

@app.callback(
    Output('choose-inv', 'children'),
    Output('inv-div', 'style'),
    Input('result', 'children'),
    State('inv-num', 'max'),
    prevent_initial_call = True
)
def set_invNum_caption(state, inv):
    return f"Choose number of inversions (max: {inv})", {'display':'inline-block', "margin-left": "20px"}

###############################################################################

@app.callback(
    Output('choose-left', 'children'),
    Output('left-div', 'style'),
    Input('result', 'children'),
    State('left-num', 'max'),
    prevent_initial_call = True
)
def set_leftNum_caption(state, left):
    return f"Choose number of left flanks (max: {left})", {'display':'inline-block', "margin-left": "20px"}

###############################################################################

@app.callback(
    Output('choose-right', 'children'),
    Output('right-div', 'style'),
    Input('result', 'children'),
    State('right-num', 'max'),
    prevent_initial_call = True
)
def set_rightNum_caption(state, right):
    return f"Choose number of right flanks (max: {right})", {'display':'inline-block', "margin-left": "20px"}

###############################################################################

@app.callback(
    Output('choose-contig', 'children'),
    Output('contig-div', 'style'),
    Input('result', 'children'),
    State('contig-num', 'max'),
    prevent_initial_call = True
)
def set_contig_caption(state, contig_max):
    return f"Choose a contig number (max: {contig_max})", {'display':'inline-block'}

###############################################################################

@app.callback(
    Output('family-dropdown', 'options'),
    Input('numbers-dropdown', 'value'),
    State('countries-dropdown', 'value'),
    prevent_initial_call = True
)
def set_family_options(specNum, spec):
    options = []
    for fam in vFam:
        temp = []
        for vir in vFam[fam]:
            if vir in specVirus[(spec, specNum)]:
                options.append(fam)
                break
    return [{'label': i, 'value': i} for i in options]

###############################################################################

@app.callback(
    Output('virus-dropdown', 'options'),
    Input('family-dropdown', 'value'),
    State('countries-dropdown', 'value'),
    State('numbers-dropdown', 'value'),
    prevent_initial_call = True
)
def set_virus_options(family, spec, specNum):
    options = []
    for v in vFam[family]:
        if v in specVirus[(spec, specNum)]:
            options.append(v)
    return [{'label': i, 'value': i} for i in options]

###############################################################################

@app.callback(
    Output('numbers-dropdown', 'options'),
    Input('countries-dropdown', 'value'), prevent_initial_call = True
)
def set_cities_options(selected_country):
    return [{'label': i, 'value': i} for i in all_regions[selected_country]]

###############################################################################

@app.callback(
    Output('virus-dropdown', 'value'),
    Input('virus-dropdown', 'options'), prevent_initial_call = True
)
def set_virus_value(available_options):
    return available_options[0]['value']

###############################################################################

@app.callback(
    Output('contig-num', 'max'),
    Output('result','children'),
    Output('matches-num', 'max'),
    Output('inv-num', 'max'),
    Output('left-num', 'max'),
    Output('right-num', 'max'),
    Input('submit', 'n_clicks'),
    Input('contig-num', 'value'),
    State('countries-dropdown', 'value'),
    State('numbers-dropdown', 'value'),
    State('virus-dropdown', 'value'),
    prevent_initial_call = True
)
def pickContig(state, cNum, spec, specNum, virus):

    cNum -= 1
    dir = label2Specimen[(spec,int(specNum))]
    pathN = '/Volumes/Data2/specimens/'+ dir + '/results/xml/' + dir +'_hits.xml'

    nodeL = []

    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(pathN, parser).getroot()
    try:
        root.xpath(f"""/root/contig[@besthit="True"]/virushit[@stitle = '{virus}']""")[0]
    except:
        return 0, "does not exist", 0,0,0,0
    contigNames = root.xpath(f"""/root/contig[@besthit="True"][virushit[@stitle = '{virus}']]/@name""")
    #for contigName in contigNames:
    contig = root.xpath(f"""/root/contig[@besthit="True" and @name = '{contigNames[cNum]}'][virushit[@stitle = '{virus}']]""")[0]
    caption = contigNames[cNum]#f"{virus} insertion sites in {spec} {specNum} ({contigName})"
    nodeL.append(caption)
    flanks = contig.xpath('./flanks')[0]
    r = flanks.xpath("./*/@rightid")
    l = flanks.xpath("./*/@leftid")
    matches = len(flanks.xpath('./match'))
    inversions = len(flanks.findall('inversion'))
    hitsLeft = contig.xpath('./vectorhitleft/@id')
    hitsRight = contig.xpath('./vectorhitright/@id')
    rCount = 0
    lCount = 0
    for x in hitsLeft:
        if x not in l:
            lCount += 1
    for y in hitsRight:
        if y not in r:
            rCount += 1
    return len(contigNames), pathN, matches, inversions, lCount, rCount

###############################################################################

@app.callback(
    Output('contigDiagram', 'src'),
    Output('left-scroll', 'disabled'),
    Output('right-scroll', 'disabled'),
    Output('buttons', 'style'),
    Output('download-seq', 'children'),
    Input('result', 'children'),
    Input('contig-num', 'value'),
    Input('matches-num', 'value'),
    Input('inv-num', 'value'),
    Input('left-num', 'value'),
    Input('right-num', 'value'),
    Input('left-scroll', 'n_clicks'),
    Input('right-scroll', 'n_clicks'),
    State('virus-dropdown', 'value'),
    State('countries-dropdown', 'value'),
    State('numbers-dropdown', 'value'),
    State('matches-num', 'max'),
    State('inv-num', 'max'),
    State('left-num', 'max'),
    State('right-num', 'max'),
    prevent_initial_call = True
)
def drawContig(path, contigNum, matchNum, invNum, leftFlanks, rightFlanks, prev, next, virus, spec, specNum, mMax, invMax, lMax, rMax):

    if path == "does not exist":
        fig = pyplot.figure()
        fig.text(.5, .5, 'No data found', ha='center')
        buf = io.BytesIO()
        pyplot.savefig(buf, format='png')
        pyplot.close()
        data = base64.b64encode(buf.getbuffer()).decode('utf8')
        pic = f"""data:image/png;base64,{data}"""
        return pic, True, True, {'display': 'block', 'textAlign': 'center'}
    chromNames = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}
    totalMatch = 0
    totalInv = 0
    totalLeft = 0
    totalRight = 0
    labelFontSize = 6
    axesFontSize = 8
    thickness = 12

    contigNum -= 1

    imageNum = next - prev
    flankNum = 0

    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(path, parser).getroot()

    contigNames = root.xpath(f"""/root/contig[@besthit="True"][virushit[@stitle = '{virus}']]/@name""")
    node = contigNames[contigNum]

    contig = root.xpath(f"""/root/contig[@name="{node}"][@besthit="True"]""")[0]

    virusHits = contig.xpath("""./virushit""")
    features = []
    newline = "\n"
    sequence =  ""
    contigLength = int(contig.xpath("""./@name""")[0].split('_')[3])
    for virusHit in virusHits:
        qstart = int(virusHit.xpath("""./qstart/text()""")[0])
        qend = int(virusHit.xpath("""./qend/text()""")[0])
        qseq = virusHit.xpath("""./qseq/text()""")[0]
        sstart = int(virusHit.xpath("""./sstart/text()""")[0])
        send = int(virusHit.xpath("""./send/text()""")[0])
        evalue = virusHit.xpath("""./evalue/text()""")[0]
        bitscore = virusHit.xpath("""./bitscore/text()""")[0]
        pident = virusHit.xpath("""./pident/text()""")[0]
        stitle = virusHit.xpath("""./@stitle""")[0]
        seqid = virusHit.xpath("""./@seqid""")[0][3:-1]
        length = abs(sstart - send) + 1
        if send > sstart:
            strandStr = '+'
        else:
            strandStr = '-'

        if qend > qstart:
            strand = 1
        else:
            strand = -1
        gf = GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                               color='#ff0000', fontdict = {'size': labelFontSize},
                                               label = stitle.lstrip('|') + ' ' + str(sstart) + '-' + str(send) + ' (' + str(length) + ' bp; ' + str(pident) + '%)')  #; ' + str(evalue) + ')')
        features.append(gf)
        sequence +=  f">{node}|{qstart}-{qend}{newline}{qseq}{newline}"

    caption = f"{virus} insertion sites in {spec} {specNum} ({node})"
    hitsLeft = contig.xpath('./vectorhitleft')
    hitsRight = contig.xpath('./vectorhitright')
    hitsOverlap = contig.xpath('./vectorhitoverlap')
    vqstart = int(contig.xpath('./virushit/qstart/text()')[0])
    vqend = int(contig.xpath('./virushit/qend/text()')[-1])
    flanks = contig.xpath('./flanks')[0]

    r = flanks.xpath("./*/@rightid")
    l = flanks.xpath("./*/@leftid")

    hitsDone = l + r

    rCount = 0
    lCount = 0
    for x in contig.xpath("./vectorhitleft/@id"):
        if x not in l:
            lCount += 1
    for y in contig.xpath("./vectorhitright/@id"):
        if y not in r:
            rCount += 1

    matches = flanks.xpath('./match')

    if matchNum > 0:
        totalMatch = math.ceil(len(matches)/matchNum)
        if len(matches) > 0:
            matchCount = matchNum*imageNum +1
            for match in matches[imageNum*matchNum : imageNum*matchNum + matchNum]:
                for hits, attribName in [(hitsLeft, 'leftid'), (hitsRight, 'rightid')]:
                    for v in hits:
                        if v.attrib['id'] == match.attrib[attribName]:
                            hitsDone.append(v.attrib['id'])
                            qstart = int(v.find('qstart').text)
                            qend = int(v.find('qend').text)
                            if qend > qstart:
                                strand = 1
                            else:
                                strand = -1
                            sstart = int(v.find('sstart').text)
                            send = int(v.find('send').text)
                            if send > sstart:
                                strandStr = '+'
                            else:
                                strandStr = '-'
                            seqid = v.attrib['seqid']
                            if seqid in chromNames:
                                seqid = chromNames[seqid]
                            elif seqid[:8] == 'NW_01873':
                                seqid = 'Scaffold ' + seqid[8:12]
                            if attribName == 'leftid':
                                distance = vqstart - qend
                            else:
                                distance = qstart - vqend
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                                       color='#009051', fontdict = {'size': labelFontSize},
                                                       label = 'M' + str(matchCount) + ': ' + seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + '), Dist: ' + str(distance) + ' bp'))
                            break
                matchCount += 1

    inversions = flanks.findall('inversion')

    if invNum > 0:
        totalInv = math.ceil(len(inversions)/invNum)
        if len(inversions) > 0:
            matchCount = invNum*imageNum +1
            for match in inversions[invNum*imageNum : invNum*imageNum+ invNum]:
                for hits, attribName in [(hitsLeft, 'leftid'), (hitsRight, 'rightid')]:
                    for v in hits:
                        if v.attrib['id'] == match.attrib[attribName]:
                            hitsDone.append(v.attrib['id'])
                            qstart = int(v.find('qstart').text)
                            qend = int(v.find('qend').text)
                            if qend > qstart:
                                strand = 1
                            else:
                                strand = -1
                            sstart = int(v.find('sstart').text)
                            send = int(v.find('send').text)
                            if send > sstart:
                                strandStr = '+'
                            else:
                                strandStr = '-'
                            seqid = v.attrib['seqid']
                            if seqid in chromNames:
                                seqid = chromNames[seqid]
                            elif seqid[:8] == 'NW_01873':
                                seqid = 'Scaffold ' + seqid[8:12]
                            if attribName == 'leftid':
                                distance = vqstart - qend
                            else:
                                distance = qstart - vqend
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                                       color='indigo', fontdict = {'size': labelFontSize},
                                                       label = 'I' + str(matchCount) + ': ' + seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + '), Dist: ' + str(distance) + ' bp'))
                            break
                matchCount += 1



    aaCoverage = [0] * contigLength
    for hits in (hitsLeft, hitsRight): #, hitsOverlap):
        if len(hits) > 0:
            hitsDict = []  # consolidate query hit: [subject hits]
            countS = 0
            for v in hits:
                pos = (v.find('qstart').text, v.find('qend').text)
                for x in range(int(pos[0]), int(pos[1])):
                    aaCoverage[x] += 1
                if v.attrib['id'] not in hitsDone:
                    hitsDict.append((v.attrib['id'], pos, [v.find('sstart').text, v.find('send').text, v.attrib['seqid']]))

            if hits != hitsRight:
                hitsDict.sort(key = lambda tup: tup[1][1], reverse = True)  # closest first
            else:
                hitsDict.sort(key = lambda tup: tup[1][1])  # closest first

            if leftFlanks > 0 or rightFlanks > 0:
                if hits == hitsLeft:
                    flankNum = leftFlanks
                    if leftFlanks != 0:
                        totalLeft = math.ceil(lCount/leftFlanks)
                elif hits == hitsRight:
                    flankNum = rightFlanks
                    if rightFlanks != 0:
                        totalRight = math.ceil(rCount/rightFlanks)

                for q in hitsDict[flankNum*imageNum : flankNum*imageNum + flankNum]:
                    qstart = int(q[1][0])
                    qend = int(q[1][1])
                    if qend > qstart:
                        strand = 1
                    else:
                        strand = -1
                    s = q[2]
                    sstart = int(s[0])
                    send = int(s[1])
                    if send > sstart:
                        strandStr = '+'
                    else:
                        strandStr = '-'
                    seqid = s[2]
                    if seqid in chromNames:
                        seqid = chromNames[seqid]
                    elif seqid[:8] == 'NW_01873':
                        seqid = 'Scaffold ' + seqid[8:12]
                    if hits == hitsOverlap:
                        color = '#00ff00'
                    else:
                        color = '#0000ff'
                    if hits == hitsLeft:
                        distance = vqstart - qend
                    else:
                        distance = qstart - vqend
                    features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                               color=color, fontdict = {'size': labelFontSize},
                                               label = seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + '), Dist: ' + str(distance) + ' bp'))
    record = GraphicRecord(sequence_length = contigLength, features = features)


    fig, (ax1, ax2) = pyplot.subplots(2, 1, sharex = True, figsize = (10, 6), gridspec_kw = {'height_ratios': [5, 1]})

    record.plot(max_label_length = 80, ax = ax1, with_ruler = False)
    ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
    ax2.tick_params(axis='both', which='major', labelsize=axesFontSize)
    ax2.set_ylim(bottom = 0, top = max(aaCoverage + [1]))
    ax2.set_yticks([0, max(aaCoverage + [1]) // 2, max(aaCoverage + [1])])
    ax2.set_ylabel('Aa hits', fontsize = axesFontSize)
    fig.text(.5, .01, caption, ha='center')

    buf = io.BytesIO()
    pyplot.savefig(buf, format='jpg')
    pyplot.close()
    data = base64.b64encode(buf.getbuffer()).decode('utf8')
    pic = f"""data:image/png;base64,{data}"""

    imgmax = max(totalMatch, totalInv, totalLeft, totalRight) - 1

    if imgmax <= 0:
        out = [True, True]
    elif imageNum <= 0:
        out = [True, False]
    elif imageNum >= imgmax:
        out = [False, True]
    else:
        out = [False, False]

    return pic, out[0], out[1], {'display': 'block', 'textAlign': 'center'}, sequence

if __name__ == '__main__':
    app.run_server(debug=True)
