import streamlit as st

from redox_reaction import RedoxReaction

st.title("Afstem en redoxreaktion")

st.markdown(
    """
    Dette program er lavet med det udgangspunkt, at H altid har
    oxidationstallet +1 og at O altid har oxidationstallet -2. Denne regel
    bliver for eksempel brudt ved peroxider, hvor O har oxidationstallet
    -1, eller i metalforbindelser, hvor H har oxidationstallet -1.
    I sådanne tilfælde tager programmet ***IKKE*** højde for
    uregelmæssigheden og vil sandsynligvis give en fejl.

    Da programmet adskiller molekylerne fra hinanden med tegnet "+",
    kan "+"-tegnet ikke også bruges til at beskrive molekylets ladning.
    Ladningerne betegnes derfor med to "selvopfundne grundstoffer": "Lp"
    for "Ladning positiv" og "Ln" for "Ladning negativ".

    Det vil sige, at fx Fe<sup>2+</sup> bliver til "FeLp2",
    mens OH<sup>-</sup> bliver til "OHLn1".
    """,
    unsafe_allow_html=True
)

unbalanced_equation = st.text_input(
    label="Indtast en redoxreaktion ifølge ovennævnte form",
    value="FeLp2 + MnO4Ln1 -> FeLp3 + MnLp2"
)

ph_choice = st.radio(
    label="Vælg om reaktionen foregår i et surt, basisk \
            eller neutralt miljø",
    options=("Surt", "Basisk", "Neutralt")
)

if ph_choice == "Surt":
    ph = "a"
elif ph_choice == "Basisk":
    ph = "b"
else:
    ph = "n"

equation = RedoxReaction(unbalanced_equation, ph)

balanced_equation = equation.balance()
unbalanced_equation_output = equation.format_unbalanced_equation()

st.header("Ikke-afstemt reaktion:")

st.markdown(
    "**" + unbalanced_equation_output + "**",
    unsafe_allow_html=True
)

st.header("Afstemt reaktion:")

st.markdown(
    "**" + balanced_equation + "**",
    unsafe_allow_html=True
)
