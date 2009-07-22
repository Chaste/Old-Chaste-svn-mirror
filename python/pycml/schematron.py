#!/usr/bin/env python
#Warning: this is an auto-generated file.  Do not edit unless you're
#sure you know what you're doing

import sys
import codecs
import optparse
import cStringIO
from Ft.Xml import InputSource, CreateInputSource
from Ft.Xml.Xslt import PatternList, parser
from Ft.Xml.Xslt import Stylesheet, Processor, OutputHandler, OutputParameters
from Ft.Xml.XPath import Compile as CompileXPath
from Ft.Xml.XPath.Context import Context as XPathContext
from Ft.Xml.Xslt.XsltContext import XsltContext
from Ft.Xml.Domlette import NonvalidatingReader, GetAllNs
from Ft.Xml.XPath import Conversions
from Ft.Xml.XPath import Util
from Ft.Xml.XPath import CoreFunctions

from Ft.Xml.Xslt.XmlWriter import XmlWriter
from Ft.Xml.Xslt.OutputParameters import OutputParameters
from Ft.Xml.Xslt import XsltFunctions, Exslt
from Ft.Lib import Uri

from amara import domtools

XPATTERN_PARSER = parser.new()
del parser

#STRON_DOC_DUMMY_BASE = 'http://4suite.org/amara/scimitar/schematron-file'

STRON_BASE_URI = 'urn:uuid:c64628c2-f75b-49bf-abf6-2f09b16504bb'
QUERY_BINDING = u'xslt'


def rule1(node):
    #For context XPattern u'cellml:component'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'not(@name=following-sibling::cellml:component/@name)')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Component names must be unique within a model (3.4.2.2).\n        The name '")
        expr = CompileXPath(u'@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is repeated.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule2(node):
    #For context XPattern u'cellml:component/cellml:variable'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'not(@name=following-sibling::cellml:variable/@name)')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Variable names must be unique within a component (3.4.3.2).\n        The name '")
        expr = CompileXPath(u'@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is repeated in component '")
        expr = CompileXPath(u'../@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule3(node):
    #For context XPattern u'cellml:units'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'u')] = CompileXPath(u'@name').evaluate(xpath_ctx)
    expr = CompileXPath(u'not(@name=following-sibling::cellml:units/@name)')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Units names must be unique within the parent component or\n        model (5.4.1.2).\n        The name '")
        expr = CompileXPath(u'@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is repeated in '")
        expr = CompileXPath(u'../@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n        (Note however that a definition in a component may override\n         a definition in the parent model.)\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u"not($u='ampere' or $u='becquerel' or $u='candela' or $u='celsius' or $u='coulomb' or $u='dimensionless' or $u='farad' or $u='gram' or $u='gray' or $u='henry' or $u='hertz' or $u='joule' or $u='katal' or $u='kelvin' or $u='kilogram' or $u='liter' or $u='litre' or $u='lumen' or $u='lux' or $u='meter' or $u='metre' or $u='mole' or $u='newton' or $u='ohm' or $u='pascal' or $u='radian' or $u='second' or $u='siemens' or $u='sievert' or $u='steradian' or $u='tesla' or $u='volt' or $u='watt' or $u='weber')")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        The standard unit '")
        expr = CompileXPath(u'$u')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' may not be redefined (5.4.1.2).\n        This is attempted in '")
        expr = CompileXPath(u'../@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule4(node):
    #For context XPattern u'cellml:variable_ref'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'not(@variable=following-sibling::cellml:variable_ref/@variable)')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        A variable must only be referenced once in a single reaction\n        (7.4.2.2).\n        The variable '")
        expr = CompileXPath(u'@variable')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is referenced repeatedly,\n        in a reaction in component '")
        expr = CompileXPath(u'ancestor::cellml:component/@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule5(node):
    #For context XPattern u'cellml:role'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'dv')] = CompileXPath(u'@delta_variable').evaluate(xpath_ctx)
    expr = CompileXPath(u'not($dv) or count(ancestor::cellml:component/cellml:reaction//cellml:role[@delta_variable=$dv]) < 2')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        The value of the delta_variable attribute must be unique across\n        all role elements contained within the parent component element\n        (7.4.3.7).\n        The variable '")
        expr = CompileXPath(u'$dv')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is repeated in component '")
        expr = CompileXPath(u'ancestor::cellml:component/@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule6(node):
    #For context XPattern u'cellml:component/cellml:variable'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'u')] = CompileXPath(u'@units').evaluate(xpath_ctx)
    expr = CompileXPath(u"$u=ancestor::*/cellml:units/@name or $u='ampere' or $u='becquerel' or $u='candela' or $u='celsius' or $u='coulomb' or $u='dimensionless' or $u='farad' or $u='gram' or $u='gray' or $u='henry' or $u='hertz' or $u='joule' or $u='katal' or $u='kelvin' or $u='kilogram' or $u='liter' or $u='litre' or $u='lumen' or $u='lux' or $u='meter' or $u='metre' or $u='mole' or $u='newton' or $u='ohm' or $u='pascal' or $u='radian' or $u='second' or $u='siemens' or $u='sievert' or $u='steradian' or $u='tesla' or $u='volt' or $u='watt' or $u='weber'")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        The value of the units attribute on a variable must be either\n        one of the standard units or the name of a unit defined in the\n        current component or model (3.4.3.3).\n        The units '")
        expr = CompileXPath(u'$u')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' on the variable '")
        expr = CompileXPath(u'@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in component '")
        expr = CompileXPath(u'../@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' do not qualify.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule7(node):
    #For context XPattern u'cellml:map_components'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'c1')] = CompileXPath(u'@component_1').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'c2')] = CompileXPath(u'@component_2').evaluate(xpath_ctx)
    expr = CompileXPath(u'$c1=/cellml:model/cellml:component/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Connections must be between components defined in the current model\n        (3.4.5.2).\n        There is no component '")
        expr = CompileXPath(u'$c1')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u'$c2=/cellml:model/cellml:component/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Connections must be between components defined in the current model\n        (3.4.5.3).\n        There is no component '")
        expr = CompileXPath(u'$c2')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u'$c1 != $c2')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        A connection must link two different components (3.4.5.4).\n        The component '")
        expr = CompileXPath(u'$c1')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is being connected to itself.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u'not(following::cellml:map_components[@component_1=$c1 and @component_2=$c2 or @component_1=$c2 and @component_2=$c1])')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Each map_components element must map a unique pair of components\n        (3.4.5.4).\n        The pair ('")
        expr = CompileXPath(u'$c1')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"','")
        expr = CompileXPath(u'$c2')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"') is repeated.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule8(node):
    #For context XPattern u'cellml:map_variables'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'c1')] = CompileXPath(u'../cellml:map_components/@component_1').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'c2')] = CompileXPath(u'../cellml:map_components/@component_2').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'v1')] = CompileXPath(u'@variable_1').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'v2')] = CompileXPath(u'@variable_2').evaluate(xpath_ctx)
    expr = CompileXPath(u'$v1=/cellml:model/cellml:component[@name=$c1]/cellml:variable/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        A variable mapping must be between existing variables (3.4.6.2).\n        Variable '")
        expr = CompileXPath(u'$v1')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' doesn't exist in component '")
        expr = CompileXPath(u'$c1')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u'$v2=/cellml:model/cellml:component[@name=$c2]/cellml:variable/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        A variable mapping must be between existing variables (3.4.6.3).\n        Variable '")
        expr = CompileXPath(u'$v2')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' doesn't exist in component '")
        expr = CompileXPath(u'$c2')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule9(node):
    #For context XPattern u'mathml:ci'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'normalize-space(text())=ancestor::cellml:component/cellml:variable/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        The content of a MathML ci element must match the name of a\n        variable in the enclosing component, once whitespace normalisation\n        has been performed (4.4.2.1).\n        Variable '")
        expr = CompileXPath(u'text()')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' does not exist in component '")
        expr = CompileXPath(u'ancestor::cellml:component/@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule10(node):
    #For context XPattern u'mathml:cn'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'u')] = CompileXPath(u'@cellml:units').evaluate(xpath_ctx)
    expr = CompileXPath(u"$u=ancestor::*/cellml:units/@name or $u='ampere' or $u='becquerel' or $u='candela' or $u='celsius' or $u='coulomb' or $u='dimensionless' or $u='farad' or $u='gram' or $u='gray' or $u='henry' or $u='hertz' or $u='joule' or $u='katal' or $u='kelvin' or $u='kilogram' or $u='liter' or $u='litre' or $u='lumen' or $u='lux' or $u='meter' or $u='metre' or $u='mole' or $u='newton' or $u='ohm' or $u='pascal' or $u='radian' or $u='second' or $u='siemens' or $u='sievert' or $u='steradian' or $u='tesla' or $u='volt' or $u='watt' or $u='weber'")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Units on a cn element must be standard or defined in the current\n        component or model (4.4.3.2).\n        Units '")
        expr = CompileXPath(u'$u')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' are not defined in component '")
        expr = CompileXPath(u'ancestor::cellml:component/@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule11(node):
    #For context XPattern u'cellml:unit'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'u')] = CompileXPath(u'@units').evaluate(xpath_ctx)
    expr = CompileXPath(u"$u=ancestor::*/cellml:units/@name or $u='ampere' or $u='becquerel' or $u='candela' or $u='celsius' or $u='coulomb' or $u='dimensionless' or $u='farad' or $u='gram' or $u='gray' or $u='henry' or $u='hertz' or $u='joule' or $u='katal' or $u='kelvin' or $u='kilogram' or $u='liter' or $u='litre' or $u='lumen' or $u='lux' or $u='meter' or $u='metre' or $u='mole' or $u='newton' or $u='ohm' or $u='pascal' or $u='radian' or $u='second' or $u='siemens' or $u='sievert' or $u='steradian' or $u='tesla' or $u='volt' or $u='watt' or $u='weber'")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        The value of the units attribute on a unit element must be taken\n        from the dictionary of standard units or be the name of a \n        user-defined unit in the current component or model (5.4.2.2).\n        Units '")
        expr = CompileXPath(u'$u')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' are not defined.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule12(node):
    #For context XPattern u'cellml:relationship_ref'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'not(@relationship=following-sibling::cellml:relationship_ref/@relationship and @name=following-sibling::cellml:relationship_ref/@name)')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n         A group element must not contain two or more relationship_ref\n         elements that define a relationship attribute in a common\n         namespace with the same value and that have the same name\n         attribute value (which may be non-existent) (6.4.2.5).\n         Relationship '")
        expr = CompileXPath(u'@relationship')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' name '")
        expr = CompileXPath(u'@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is repeated.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule13(node):
    #For context XPattern u"cellml:group[cellml:relationship_ref/@relationship='encapsulation' or cellml:relationship_ref/@relationship='containment']/cellml:component_ref"
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'cellml:component_ref')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Containment and encapsulation relationships must be hierarchical\n        (6.4.3.2).\n        Potentially top-level component '")
        expr = CompileXPath(u'@component')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'\n        has not been given children.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule14(node):
    #For context XPattern u'cellml:component_ref'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'hierarchies')] = CompileXPath(u'ancestor::cellml:group/cellml:relationship_ref/@relationship').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'c')] = CompileXPath(u'@component').evaluate(xpath_ctx)
    expr = CompileXPath(u'not(cellml:component_ref) or not(following::cellml:component_ref[@component=$c and $hierarchies=ancestor::cellml:group/cellml:relationship_ref/@relationship]/cellml:component_ref)')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        In a given hierarchy, only one component_ref element referencing\n        a given component may contain children (6.4.3.2).\n        Component '")
        expr = CompileXPath(u'@component')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' has children in multiple locations.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u'not(ancestor::cellml:component_ref) or not(following::cellml:component_ref[$hierarchies=ancestor::cellml:group/cellml:relationship_ref/@relationship]/cellml:component_ref[@component=$c])')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        In a given hierarchy, a component may not be a child more than\n        once (6.4.3.2).\n        Component '")
        expr = CompileXPath(u'@component')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' has multiple parents.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u'$c=/cellml:model/cellml:component/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        A component_ref element must reference a component in the current\n        model (6.4.3.3).\n        Component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' does not exist.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule15(node):
    #For context XPattern u'cellml:role'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    xpath_ctx.varBindings[(None, u'c')] = CompileXPath(u'ancestor::cellml:component/@name').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'v')] = CompileXPath(u'../@variable').evaluate(xpath_ctx)
    expr = CompileXPath(u"not(/cellml:model/cellml:component_ref[@component=$c]/cellml:component_ref[ancestor::cellml:group/cellml:relationship_ref/@relationship='encapsulation']) or not(@delta_variable)")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        An encapsulating component must not contain delta_variable\n        attrs on any role elements (7.4.1.3).\n        Variable '")
        expr = CompileXPath(u'$v')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' has a role\n        which specifies a delta_variable.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u"not(/cellml:model/cellml:component_ref[@component=$c]/cellml:component_ref[ancestor::cellml:group/cellml:relationship_ref/@relationship='encapsulation']) or not(mathml:math)")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        An encapsulating component must not contain explicit mathematics\n        within role elements (7.4.1.3).\n        Variable '")
        expr = CompileXPath(u'$v')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' has a role\n        with explicit mathematics.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u"not(@role='rate' and count(ancestor::cellml:reaction//cellml:role[@role='rate']) > 1)")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        There may only be one rate variable per reaction (7.4.3.3).\n        Variable '")
        expr = CompileXPath(u'$v')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in a reaction in component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'\n        is one of multiple rate variables.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u"not(ancestor::cellml:reaction/@reversible='no') or (@direction='forward' or not(@direction))")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        Only reversible reactions may occur in 2 directions (7.4.3.5).\n        Variable '")
        expr = CompileXPath(u'$v')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in a non-reversible reaction in\n        component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' specifies a non-forward direction.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    xpath_ctx.varBindings[(None, u'role')] = CompileXPath(u'@role').evaluate(xpath_ctx)
    xpath_ctx.varBindings[(None, u'dir')] = CompileXPath(u'@direction').evaluate(xpath_ctx)
    expr = CompileXPath(u'not(following-sibling::cellml:role[@role=$role and @direction=$dir])')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n         Each role element contained in a given variable_ref element\n         must have a unique combination of values for the role and\n         direction attributes (7.4.3.5).\n         Variable '")
        expr = CompileXPath(u'$v')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' has conflicting roles.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')
    expr = CompileXPath(u"not(@delta_variable and @stoichiometry) or ../../cellml:variable_ref/cellml:role[@role='rate']")
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n         If the delta_variable and stoichiometry attributes are both\n         declared on any single reaction participant, a variable must\n         be provided to represent the reaction rate (7.4.3.8).\n         The variable '")
        expr = CompileXPath(u'$v')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' in component '")
        expr = CompileXPath(u'$c')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' has\n         a role with delta_variable and stoichiometry set, but there is\n         no rate variable.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


def rule16(node):
    #For context XPattern u'variable_ref'
    vars = {}
    xpath_ctx = XsltContext(node, processor=key_handler,
                            processorNss=NSS, varBindings=vars,
                            currentNode=node
                            )
    xpath_ctx.currentInstruction = faux_instruction()
    xpath_ctx.functions = FUNCTIONS
    expr = CompileXPath(u'@variable=ancestor::cellml:component/cellml:variable/@name')
    if not Conversions.BooleanValue(expr.evaluate(xpath_ctx)):
        WRITER.text(u'Assertion failure:\n')
        WRITER.text(u"\n        A variable reference in a reaction must reference a variable\n        defined in the current component (7.4.2.2).\n        Variable '")
        expr = CompileXPath(u'@variable')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"' is not defined in '")
        expr = CompileXPath(u'ancestor::cellml:component/@name')
        val = Conversions.StringValue(expr.evaluate(xpath_ctx))
        WRITER.text(val)
        WRITER.text(u"'.\n      ")
        WRITER.text(u'\n')
        WRITER.text(u'\n')

    return


PATTERNS = [
(
 u'Checking uniqueness of names', None, {
  XPATTERN_PARSER.parse(u'cellml:component'): rule1,
  XPATTERN_PARSER.parse(u'cellml:units'): rule3,
  XPATTERN_PARSER.parse(u'cellml:component/cellml:variable'): rule2,
  XPATTERN_PARSER.parse(u'cellml:variable_ref'): rule4,
  XPATTERN_PARSER.parse(u'cellml:role'): rule5,
}),
(
 u'Rules from section 3', None, {
  XPATTERN_PARSER.parse(u'cellml:map_variables'): rule8,
  XPATTERN_PARSER.parse(u'cellml:component/cellml:variable'): rule6,
  XPATTERN_PARSER.parse(u'cellml:map_components'): rule7,
}),
(
 u'Rules from section 4', None, {
  XPATTERN_PARSER.parse(u'mathml:cn'): rule10,
  XPATTERN_PARSER.parse(u'mathml:ci'): rule9,
}),
(
 u'Rules from section 5', None, {
  XPATTERN_PARSER.parse(u'cellml:unit'): rule11,
}),
(
 u'Rules from section 6', None, {
  XPATTERN_PARSER.parse(u'cellml:component_ref'): rule14,
  XPATTERN_PARSER.parse(u'cellml:relationship_ref'): rule12,
  XPATTERN_PARSER.parse(u"cellml:group[cellml:relationship_ref/@relationship='encapsulation' or cellml:relationship_ref/@relationship='containment']/cellml:component_ref"): rule13,
}),
(
 u'Rules from section 7', None, {
  XPATTERN_PARSER.parse(u'cellml:role'): rule15,
  XPATTERN_PARSER.parse(u'variable_ref'): rule16,
}),
]

PHASES = {
}

CONTEXTS = {
}

DIAGNOSTICS = {
}

NSS = {
u'mathml': u'http://www.w3.org/1998/Math/MathML', u'cmeta': u'http://www.cellml.org/metadata/1.0#', u'dc': u'http://purl.org/dc/elements/1.1/', u'rdf': u'http://www.w3.org/1999/02/22-rdf-syntax-ns#', u'bqs': u'http://www.cellml.org/bqs/1.0#', u'cellml': u'http://www.cellml.org/cellml/1.0#', u'dcterms': u'http://purl.org/dc/terms/', u'vCard': u'http://www.w3.org/2001/vcard-rdf/3.0#', }

KEYS = [
]

DEFAULT_PHASE = u'#ALL'


#Set up the function library for context objects according to queryBinding attr
#Determine from the query binding whether functions could require a result tree
RESULT_TREE_PROCESSOR = False
if QUERY_BINDING == u'full':
    FUNCTIONS = XsltContext.functions.copy()
    RESULT_TREE_PROCESSOR = True
elif QUERY_BINDING == u'exslt':
    FUNCTIONS = CoreFunctions.CoreFunctions.copy()
    FUNCTIONS.update(XsltFunctions.CoreFunctions)
    FUNCTIONS.update(Exslt.ExtFunctions)
    RESULT_TREE_PROCESSOR = True
elif QUERY_BINDING == u'xslt':
    FUNCTIONS = CoreFunctions.CoreFunctions.copy()
    FUNCTIONS.update(XsltFunctions.CoreFunctions)
elif QUERY_BINDING == u'xpath':
    FUNCTIONS = CoreFunctions.CoreFunctions.copy()
else:
    raise ValueError('Unknown query binding: %s'%repr(QUERY_BINDING))


class faux_instruction:
    baseUri = STRON_BASE_URI

class faux_root_node:
    sources = []

class faux_xslt_proc(Stylesheet.StylesheetElement, Processor.Processor):
    #Pretends to be an XSLT processor for processing keys
    def __init__(self, doc):
        self.namespaces = NSS
        #self._keys = [ (Util.ExpandQName(k[0], doc), k[1], k[2]) for k in KEYS ]
        keys = [ (Util.ExpandQName(k[0], doc), (k[1], k[2], GetAllNs(doc))) for k in KEYS ]
        self.keys = {}
        self._keys = {}
        for name, info in keys:
            self._keys.setdefault(name, []).append(info)
        #Update all the keys for all documents in the context
        #Note: 4Suite uses lazy key eval.  Consider emulating this
        for key in self._keys:
            Stylesheet.StylesheetElement.updateKey(self, doc, key[0], self)
        if RESULT_TREE_PROCESSOR:
            Processor.Processor.__init__(self)
            #Start code from Processor.runNode
            if hasattr(doc, 'baseURI'):
                node_baseUri = doc.baseURI
            elif hasattr(doc, 'refUri'):
                node_baseUri = doc.refUri
            else:
                node_baseUri = None
            sourceUri = node_baseUri or Uri.BASIC_RESOLVER.generate()
            #Create a dummy iSrc
            docInputSource = InputSource.InputSource(
                None, sourceUri, processIncludes=1,
                stripElements=self.getStripElements(),
                factory=self.inputSourceFactory)
            #return self.execute(node, docInputSource, outputStream=None)
            outputStream = cStringIO.StringIO()
            def writer_changed_handler(newWriter):
                self.writers[-1] = newWriter
            self.outputParams = OutputParameters()
            writer = OutputHandler.OutputHandler(
                self.outputParams, outputStream, writer_changed_handler)
            self.writers = [writer]
        self.initialFunctions = {}
        self.stylesheet = self
        self.root = faux_root_node()
        self._documentInputSource = InputSource.DefaultFactory.fromString('<dummy/>', Uri.OsPathToUri('.'))
        return


def validate(xmlf, reportf, phase=None):
    global WRITER, EXTFUNCTS
    oparams = OutputParameters()
    oparams.indent = 'yes'
    WRITER = XmlWriter(oparams, reportf)
    WRITER.startDocument()

    WRITER.text(u'Processing schema: ')
    WRITER.text(u'A Schematron Schema for CellML 1.0')
    WRITER.text(u'\n\n')
    
    if xmlf == '-':
        doc = NonvalidatingReader.parseStream(sys.stdin, 'urn:stron-candidate-dummy')
    elif not isinstance(xmlf, str):
        doc = NonvalidatingReader.parseStream(xmlf, 'urn:stron-candidate-dummy')
    else:
        try:
            doc = NonvalidatingReader.parseUri(xmlf)
        except ValueError:
            doc = NonvalidatingReader.parseUri(Uri.OsPathToUri(xmlf))

    global key_handler
    key_handler = faux_xslt_proc(doc)
    #Pre-process keys, if any
    #for kname, (kuse, kmatch) in KEYS.items:

    if phase in [None, u'#DEFAULT']:
        phase = DEFAULT_PHASE
    if phase == u'#ALL':
        patterns = PATTERNS
    else:
        patterns = [ p for p in PATTERNS if p[1] in PHASES[phase] ]

    #Main rule-processing loop
    if patterns:
        for pat in patterns:
            #if not pat[2].keys():
            #    WRITER.text(_(u'Pattern context not given'))
            plist = PatternList(pat[2].keys(), NSS)
            WRITER.text(u'Processing pattern: ')
            WRITER.text(pat[0] or u'[unnamed]')
            WRITER.text(u'\n\n')
            for node in domtools.doc_order_iter(doc):
                matches = plist.lookup(node)
                if matches:
                    func = pat[2][matches[0]]
                    func(node)
    else:
        #Second parameter is a dictionary of prefix to namespace mappings
        plist = PatternList(CONTEXTS.keys(), {})
        for node in domtools.doc_order_iter(doc):
            #FIXME: this is a 4Suite bug work-around.  At some point remove
            #The redundant context setting
            matches = plist.lookup(node)
            if matches:
                func = CONTEXTS[matches[0]]
                func(node)
    WRITER.endDocument()
    return


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def command_line_prep():
    from optparse import OptionParser
    usage = "%prog [options] xml-file\nxml-file is the XML file to be validated"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--val-output",
                      action="store", type="string", dest="val_output",
                      help="generate the validation report file FILE (XML external parsed entity format)", metavar="FILE")
    parser.add_option("-p", "--phase",
                      action="store", type="string",
                      help="Execute the schema in the specified phase", metavar="XML_ID")
    return parser

def main(argv=[__name__]):
    #Ideas borrowed from http://www.artima.com/forums/flat.jsp?forum=106&thread=4829
    if argv is None:
        argv = sys.argv
    try:
        try:
            optparser = command_line_prep()
            global OPTIONS, ARGS
            (OPTIONS, ARGS) = optparser.parse_args(argv)
            candidate_file = ARGS[1]
        except KeyboardInterrupt:
            pass
        except:
             raise Usage(optparser.format_help())
        enc, dec, inwrap, outwrap = codecs.lookup('utf-8')
        fout = OPTIONS.val_output
        phase = OPTIONS.phase
        if fout:
            fout = open(fout, 'w')
        else:
            fout = sys.stdout
        validate(candidate_file, outwrap(fout), phase=phase)
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))

