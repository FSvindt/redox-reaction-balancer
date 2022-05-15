"""
AFSTEM EN REDOXREAKTION

1) Tildel oxidationstal
    - H -> 1
    - O -> -2
    - grundstoffer -> 0
    - forbindelsers samlede OT -> ladningen

2) Find ændring i oxidationstal
    - skal findes for de grundstoffer, der oxideres eller reduceres

3) Afstem oxidationstal med koefficienter
    - mindste fælles multiplum

4) Afstem ladning (Lp og Ln angiver henholdsvis positiv og negativ ladning)
    - afstemmes med H+ hvis sur opløsning
    - afstemmes med OH- hvis basisk opløsning
    - afstemmes med H+ eller OH- på produktsiden hvis neutral opløsning

5) Afstem O og H med H2O
"""

import chemparse as cp
import sympy as smp

x = smp.symbols("x", real=True)


# > A class that represents a redox reaction.
class RedoxReaction:
    balanced_coefficients = {}

    def __init__(self, unbalanced_equation, ph):
        """
        The function __init__() is a special function in Python classes.
        It is known as a constructor in object oriented concepts.
        This method is called when an object is created from a class
        and it allows the class to initialize the attributes of the class.

        :param unbalanced_equation: The unbalanced equation you want to balance
        :param ph: The pH of the solution ("a", "b" or "n")
        """
        self.unbalanced_equation = unbalanced_equation
        self.ph = ph  # "a" for acid, "b" for base, "n" for neutral

    def _parse(self):
        """
        It takes a chemical equation, splits it into reactants and products,
        splits the reactants and products into individual compounds,
        and then parses the compounds into a dictionary of elements
        and their coefficients

        :return: A tuple of the reactants, products, unique elements, and all
        compounds.
        """
        reactants, products = self.unbalanced_equation.split("->")

        reactant_compounds = reactants.split("+")
        product_compounds = products.split("+")

        reactant_compounds = [c.strip() for c in reactant_compounds]
        product_compounds = [c.strip() for c in product_compounds]

        all_compounds = []
        unique_elements = {}

        for compound in reactant_compounds:
            num_elements = cp.parse_formula(compound)

            for key in num_elements:
                if key != "Lp" and key != "Ln":
                    unique_elements[key] = ""

            all_compounds.append(num_elements)

        for compound in product_compounds:
            num_elements = cp.parse_formula(compound)

            num_elements = {key: -val for key, val in num_elements.items()}

            all_compounds.append(num_elements)

        return (
            reactant_compounds,
            product_compounds,
            unique_elements,
            all_compounds
        )

        all_compounds_new = []
        for compound in all_compounds:
            compound_new = {}
            for key, value in compound.items():
                if key != "Lp" and key != "Ln":
                    compound_new[key] = value

            all_compounds_new.append(compound_new)

        return all_compounds_new

    def _get_charges(
        self,
        reactant_compounds,
        product_compounds,
        all_compounds
    ):
        """
        This function takes in a list of reactant compounds,
        a list of product compounds, and a list of all compounds,
        and returns a dictionary of charges for each compound

        :param reactant_compounds: a list of the reactant compounds
        :param product_compounds: a list of the products of the reaction
        :param all_compounds: a list of all compounds in the reaction
        :return charges: a dictionary of charges for each compound
        """
        charges = {}
        for compound_index, compound in enumerate(all_compounds):
            for key, value in compound.items():
                if value > 0:  # reactants
                    if key == "Lp":
                        charges[reactant_compounds[compound_index]] = value
                    elif key == "Ln":
                        charges[reactant_compounds[compound_index]] = -value
                    else:
                        charges[reactant_compounds[compound_index]] = 0

                if value < 0:  # products
                    if key == "Lp":
                        charges[
                            product_compounds[
                                compound_index - len(reactant_compounds)
                            ]
                        ] = value
                    elif key == "Ln":
                        charges[
                            product_compounds[
                                compound_index - len(reactant_compounds)
                            ]
                        ] = -value
                    else:
                        charges[
                            product_compounds[
                                compound_index - len(reactant_compounds)
                            ]
                        ] = 0

        return charges

    def _assign_oxidation_numbers(self):
        """
        The function assigns oxidation numbers to the
        reactants and products of a chemical reaction

        :return: A tuple of two dictionaries containing oxidation numbers
        for reactants and products
        """
        (reactant_compounds,
         product_compounds,
         _,
         all_compounds) = self._parse()

        charges = self._get_charges(
            reactant_compounds,
            product_compounds,
            all_compounds
        )

        reactant_oxidation_numbers = {}
        product_oxidation_numbers = {}

        for compound_index, compound in enumerate(all_compounds):
            reactants_on_dict = {}
            products_on_dict = {}

            # assign oxidation numbers for H and O
            for key, value in compound.items():
                if value > 0:  # reactants
                    if reactant_compounds[compound_index] in reactants_on_dict:
                        if key == "O":
                            reactants_on_dict[
                                reactant_compounds[compound_index]
                            ] += -2 * value

                        elif key == "H":
                            reactants_on_dict[
                                reactant_compounds[compound_index]
                            ] += 1 * value
                    else:
                        if key == "O":
                            reactants_on_dict[
                                reactant_compounds[compound_index]
                            ] = -2 * value

                        elif key == "H":
                            reactants_on_dict[
                                reactant_compounds[compound_index]
                            ] = 1 * value

                elif value < 0:  # products
                    if product_compounds[
                        compound_index - len(reactant_compounds)
                    ] in products_on_dict:
                        if key == "O":
                            products_on_dict[
                                product_compounds[
                                    compound_index - len(reactant_compounds)
                                ]
                            ] += -2 * abs(value)

                        elif key == "H":
                            products_on_dict[
                                product_compounds[
                                    compound_index - len(reactant_compounds)
                                ]
                            ] += 1 * abs(value)
                    else:
                        if key == "O":
                            products_on_dict[
                                product_compounds[
                                    compound_index - len(reactant_compounds)
                                ]
                            ] = -2 * abs(value)

                        elif key == "H":
                            products_on_dict[
                                product_compounds[
                                    compound_index - len(reactant_compounds)
                                ]
                            ] = 1 * abs(value)

            # assign oxidation numbers for the remaining elements
            for key, value in compound.items():
                if value > 0:  # reactants
                    compound_charge = charges[
                        reactant_compounds[compound_index]
                    ]

                    if key != "O" and key != "H" \
                            and key != "Lp" and key != "Ln":
                        if reactant_compounds[compound_index] \
                                in reactants_on_dict:
                            prev_on = reactants_on_dict[
                                reactant_compounds[compound_index]
                            ]

                            on = smp.solve(x + prev_on - compound_charge)[0]

                            reactant_oxidation_numbers[key] = on

                        else:
                            on = smp.solve(x - compound_charge)[0]

                            reactant_oxidation_numbers[key] = on

                elif value < 0:  # products
                    compound_charge = -charges[
                        product_compounds[
                            compound_index - len(reactant_compounds)
                        ]
                    ]

                    if key != "O" and key != "H" \
                            and key != "Lp" and key != "Ln":
                        if product_compounds[
                            compound_index - len(reactant_compounds)
                        ] in products_on_dict:
                            prev_on = products_on_dict[
                                product_compounds[
                                    compound_index - len(reactant_compounds)
                                ]
                            ]

                            on = smp.solve(x + prev_on - compound_charge)[0]

                            product_oxidation_numbers[key] = on

                        else:
                            on = smp.solve(x - compound_charge)[0]

                            product_oxidation_numbers[key] = on

        return (reactant_oxidation_numbers, product_oxidation_numbers)

    def _balance_oxidation_numbers(self):
        """
        The function takes the oxidation numbers of the reactants and products,
        and then finds the difference between the oxidation numbers of each
        element in the reactants and products.

        The function then finds the product of all of the differences,
        and then divides the product by each difference to find the
        coefficient that needs to be multiplied by each compound to balance
        the oxidation numbers.

        The function then multiplies each compound by the coefficient that was
        found, and then returns the balanced compounds.

        The function also updates the `balanced_coefficients` dictionary with
        the coefficients that were found.

        The function also returns the reactant compounds, product compounds,
        unique elements, and all compounds.

        :return: The return statement is returning the reactant_compounds,
        product_compounds, unique_elements, and all_compounds.
        """
        (reactant_compounds,
         product_compounds,
         unique_elements,
         all_compounds) = self._parse()

        (reactant_oxidation_numbers,
         product_oxidation_numbers) = self._assign_oxidation_numbers()

        on_differences = {}
        for key, reactant_on in reactant_oxidation_numbers.items():
            product_on = product_oxidation_numbers[key]

            on_difference = reactant_on - product_on

            on_differences[key] = on_difference

        multiple = abs(smp.prod(on_differences.values()))

        coefficients = {}
        for key, value in on_differences.items():
            coefficients[key] = abs(int(multiple) // int(value))

        for element, coeff in coefficients.items():
            for compound in all_compounds:
                if element in compound:
                    compound.update(
                        (k, v * coeff) for k, v in compound.items()
                    )

            for compound in reactant_compounds:
                if element in cp.parse_formula(compound):
                    self.balanced_coefficients[compound] = coeff

            for compound in product_compounds:
                if element in cp.parse_formula(compound):
                    self.balanced_coefficients[compound] = coeff

        return (
            reactant_compounds,
            product_compounds,
            unique_elements,
            all_compounds
        )

    def _balance_charge_if_acid(
        self,
        summed_charges,
        reactant_compounds,
        product_compounds
    ):
        """
        If the reactants have a lower charge than the products,
        add a proton to the reactants. If the products have a lower charge
        than the reactants, add a proton to the products

        :param summed_charges: A tuple of the summed charges of the
        reactants and products
        :param reactant_compounds: A list of reactant compounds
        :param product_compounds: A list of product compounds
        :return: The reactant_compounds and product_compounds lists are being
        returned.
        """
        reactants_charge, products_charge = summed_charges

        if reactants_charge < products_charge:
            reactant_compounds.append("HLp1")
            self.balanced_coefficients["HLp1"] = \
                int(products_charge - reactants_charge)
        elif reactants_charge > products_charge:
            product_compounds.append("HLp1")
            self.balanced_coefficients["HLp1"] = \
                int(reactants_charge - products_charge)

        return (reactant_compounds, product_compounds)

    def _balance_charge_if_base(
        self,
        summed_charges,
        reactant_compounds,
        product_compounds
    ):
        """
        If the reactants have a lower charge than the products,
        add a hydroxide ion to the products. If the reactants have a higher
        charge than the products, add a hydroxide ion to the reactants

        :param summed_charges: A tuple of the summed charges of the
        reactants and products
        :param reactant_compounds: A list of reactant compounds
        :param product_compounds: A list of product compounds
        :return: The reactant_compounds and product_compounds lists are being
        returned.
        """
        reactants_charge, products_charge = summed_charges

        if reactants_charge < products_charge:
            product_compounds.append("OHLn1")
            self.balanced_coefficients["OHLn1"] = \
                int(products_charge - reactants_charge)
        elif reactants_charge > products_charge:
            reactant_compounds.append("OHLn1")
            self.balanced_coefficients["OHLn1"] = \
                int(reactants_charge - products_charge)

        return (reactant_compounds, product_compounds)

    def _balance_charge_if_neutral(
        self,
        summed_charges,
        reactant_compounds,
        product_compounds
    ):
        """
        If the reactants have a net charge of less than the products, add a
        hydroxide ion to the products. If the reactants have a net charge of
        more than the products, add a proton to the products

        :param summed_charges: A tuple of the summed charges of the
        reactants and products
        :param reactant_compounds: A list of the reactant compounds
        :param product_compounds: A list of the product compounds
        :return: The reactant_compounds and product_compounds lists are being
        returned.
        """
        reactants_charge, products_charge = summed_charges

        if reactants_charge < products_charge:
            product_compounds.append("OHLn1")
            self.balanced_coefficients["OHLn1"] = \
                int(products_charge - reactants_charge)
        elif reactants_charge > products_charge:
            product_compounds.append("HLp1")
            self.balanced_coefficients["HLp1"] = \
                int(reactants_charge - products_charge)

        return (reactant_compounds, product_compounds)

    def _balance_charge(self):
        """
        Balance the equation based on whether the solution is an acid, a base,
        or neutral

        :return: The reactant_compounds, product_compounds, unique_elements,
        and all_compounds are being returned.
        """
        (reactant_compounds,
         product_compounds,
         unique_elements,
         all_compounds) = self._balance_oxidation_numbers()

        charges = self._get_charges(
            reactant_compounds,
            product_compounds,
            all_compounds
        )

        reactants_charge = 0
        products_charge = 0
        for compound in reactant_compounds:
            if compound in charges:
                reactants_charge += charges[compound]

        for compound in product_compounds:
            if compound in charges:
                products_charge += -charges[compound]

        summed_charges = (reactants_charge, products_charge)

        if self.ph == "a":  # acid
            result = self._balance_charge_if_acid(
                summed_charges,
                reactant_compounds,
                product_compounds
            )
        elif self.ph == "b":  # base
            result = self._balance_charge_if_base(
                summed_charges,
                reactant_compounds,
                product_compounds
            )
        elif self.ph == "n":  # neutral
            result = self._balance_charge_if_neutral(
                summed_charges,
                reactant_compounds,
                product_compounds
            )

        reactant_compounds, product_compounds = result

        return (
            reactant_compounds,
            product_compounds,
            unique_elements,
            all_compounds
        )

    def _balance_water(self):
        """
        If there is a net oxygen count greater than 0, add water to the
        products. If there is a net oxygen count less than 0, add water to the
        reactants

        :return: the reactant_compounds, product_compounds, unique_elements,
        and all_compounds are being returned.
        """
        (reactant_compounds,
         product_compounds,
         unique_elements,
         all_compounds) = self._balance_charge()

        oxygen_count = 0
        for compound in all_compounds:
            for element, count in compound.items():
                if element == "O":
                    oxygen_count += count

        if oxygen_count > 0:
            product_compounds.append("H2O")
            self.balanced_coefficients["H2O"] = int(abs(oxygen_count))
        elif oxygen_count < 0:
            reactant_compounds.append("H2O")
            self.balanced_coefficients["H2O"] = int(abs(oxygen_count))

        return (
            reactant_compounds,
            product_compounds,
            unique_elements,
            all_compounds
        )

    def balance(self):
        """
        The function takes the reactants and products of the reaction, and then
        creates a list of strings that represent the balanced equation.

        :return: A string of the balanced equation.
        """
        (reactant_compounds,
         product_compounds,
         _,
         _) = self._balance_water()

        balanced_equation_list = []

        for reactant in reactant_compounds:
            if reactant in self.balanced_coefficients \
                    and self.balanced_coefficients[reactant] != 1:
                balanced_equation_list.append(
                    str(self.balanced_coefficients[reactant]) + " "
                )

            split_reactant = cp.parse_formula(reactant)
            for element, count in split_reactant.items():
                if element == "Lp":
                    if count != 1:
                        balanced_equation_list.append(
                            f"<sup>{int(count)}+</sup>"
                        )
                    else:
                        balanced_equation_list.append("<sup>+</sup>")

                elif element == "Ln":
                    if count != 1:
                        balanced_equation_list.append(
                            f"<sup>{int(count)}-</sup>"
                        )
                    else:
                        balanced_equation_list.append("<sup>-</sup>")

                else:
                    balanced_equation_list.append(element)

                    if count != 1:
                        balanced_equation_list.append(
                            f"<sub>{int(count)}</sub>"
                        )

            balanced_equation_list.append(" + ")

        balanced_equation_list.pop()

        balanced_equation_list.append(" -> ")

        for product in product_compounds:
            if product in self.balanced_coefficients \
                    and self.balanced_coefficients[product] != 1:
                balanced_equation_list.append(
                    str(self.balanced_coefficients[product]) + " "
                )

            split_product = cp.parse_formula(product)
            for element, count in split_product.items():
                if element == "Lp":
                    if count != 1:
                        balanced_equation_list.append(
                            f"<sup>{int(count)}+</sup>"
                        )
                    else:
                        balanced_equation_list.append("<sup>+</sup>")

                elif element == "Ln":
                    if count != 1:
                        balanced_equation_list.append(
                            f"<sup>{int(count)}-</sup>"
                        )
                    else:
                        balanced_equation_list.append("<sup>-</sup>")

                else:
                    balanced_equation_list.append(element)

                    if count != 1:
                        balanced_equation_list.append(
                            f"<sub>{int(count)}</sub>"
                        )

            balanced_equation_list.append(" + ")

        balanced_equation_list.pop()

        balanced_equation = "".join(balanced_equation_list)

        return balanced_equation

    def format_unbalanced_equation(self):
        (reactant_compounds,
         product_compounds,
         _,
         all_compounds) = self._parse()

        unbalanced_equation_list = []

        for reactant in reactant_compounds:
            split_reactant = cp.parse_formula(reactant)
            for element, count in split_reactant.items():
                if element == "Lp":
                    if count != 1:
                        unbalanced_equation_list.append(
                            f"<sup>{int(count)}+</sup>"
                        )
                    else:
                        unbalanced_equation_list.append("<sup>+</sup>")

                elif element == "Ln":
                    if count != 1:
                        unbalanced_equation_list.append(
                            f"<sup>{int(count)}-</sup>"
                        )
                    else:
                        unbalanced_equation_list.append("<sup>-</sup>")

                else:
                    unbalanced_equation_list.append(element)

                    if count != 1:
                        unbalanced_equation_list.append(
                            f"<sub>{int(count)}</sub>"
                        )

            unbalanced_equation_list.append(" + ")

        unbalanced_equation_list.pop()

        unbalanced_equation_list.append(" -> ")

        for product in product_compounds:
            split_product = cp.parse_formula(product)
            for element, count in split_product.items():
                if element == "Lp":
                    if count != 1:
                        unbalanced_equation_list.append(
                            f"<sup>{int(count)}+</sup>"
                        )
                    else:
                        unbalanced_equation_list.append("<sup>+</sup>")

                elif element == "Ln":
                    if count != 1:
                        unbalanced_equation_list.append(
                            f"<sup>{int(count)}-</sup>"
                        )
                    else:
                        unbalanced_equation_list.append("<sup>-</sup>")

                else:
                    unbalanced_equation_list.append(element)

                    if count != 1:
                        unbalanced_equation_list.append(
                            f"<sub>{int(count)}</sub>"
                        )

            unbalanced_equation_list.append(" + ")

        unbalanced_equation_list.pop()

        unbalanced_equation = "".join(unbalanced_equation_list)

        return unbalanced_equation

# ! issue: all_compounds is a list of dicts - not a dict!
