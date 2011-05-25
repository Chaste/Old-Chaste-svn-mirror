/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "XmlTools.hpp"

#include <iostream>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/QName.hpp>
#include <xercesc/util/XMLUniDefs.hpp> // chLatin_*
#include <xercesc/framework/Wrapper4InputSource.hpp>
#include <xercesc/validators/common/Grammar.hpp>

#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/xml/dom/bits/error-handler-proxy.hxx>
#include <xsd/cxx/tree/exceptions.hxx>

#include "Exception.hpp"



xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> XmlTools::ReadXmlFile(
    const std::string& rFileName,
    const ::xsd::cxx::tree::properties<char>& rProps)
{
    xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> p_doc;
    try
    {
        // Initialise Xerces
        xercesc::XMLPlatformUtils::Initialize();
        // Set up an error handler
        ::xsd::cxx::tree::error_handler<char> error_handler;
        // Parse XML to DOM
        p_doc = XmlTools::ReadFileToDomDocument(rFileName, error_handler, rProps);
        // Any errors?
        error_handler.throw_if_failed< ::xsd::cxx::tree::parsing<char> >();
    }
    catch (const ::xsd::cxx::tree::parsing<char>& e)
    {
        Finalize();
        // Test for missing schema/xml file
#if (XSD_INT_VERSION >= 3000000L)
        const ::xsd::cxx::tree::diagnostics<char>& diags = e.diagnostics();
        const ::xsd::cxx::tree::error<char>& first_error = diags[0];
#else
        const ::xsd::cxx::tree::errors<char>& errors = e.errors();
        const ::xsd::cxx::tree::error<char>& first_error = errors[0];
#endif
        if (first_error.line() == 0u)
        {
            std::cerr << first_error << std::endl;
            EXCEPTION("Missing file parsing configuration file: " + rFileName);
        }
        else
        {
            std::cerr << e << std::endl;
            EXCEPTION("XML parsing error in configuration file: " + rFileName);
        }
    }
    catch (...)
    {
        Finalize();
        throw;
    }
    return p_doc;
}


void XmlTools::Finalize()
{
    xercesc::XMLPlatformUtils::Terminate();
}


xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> XmlTools::ReadFileToDomDocument(
        const std::string& rFileName,
        ::xsd::cxx::xml::error_handler<char>& rErrorHandler,
        const ::xsd::cxx::tree::properties<char>& rProps)
{
    using namespace xercesc;
    namespace xml = xsd::cxx::xml;

    // Get an implementation of the Load-Store (LS) interface.
    const XMLCh ls_id [] = {chLatin_L, chLatin_S, chNull};
    DOMImplementation* p_impl(DOMImplementationRegistry::getDOMImplementation(ls_id));

#if _XERCES_VERSION >= 30000
    // Xerces-C++ 3.0.0 and later.
    xml::dom::auto_ptr<DOMLSParser> p_parser(p_impl->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0));
    DOMConfiguration* p_conf(p_parser->getDomConfig());

    // Discard comment nodes in the document.
    p_conf->setParameter(XMLUni::fgDOMComments, false);

    // Enable datatype normalization.
    p_conf->setParameter(XMLUni::fgDOMDatatypeNormalization, true);

    // Do not create EntityReference nodes in the DOM tree.  No
    // EntityReference nodes will be created, only the nodes
    // corresponding to their fully expanded substitution text
    // will be created.
    p_conf->setParameter(XMLUni::fgDOMEntities, false);

    // Perform namespace processing.
    p_conf->setParameter(XMLUni::fgDOMNamespaces, true);

    // Do not include ignorable whitespace in the DOM tree.
    p_conf->setParameter(XMLUni::fgDOMElementContentWhitespace, false);

    // Enable validation.
    p_conf->setParameter(XMLUni::fgDOMValidate, true);
    p_conf->setParameter(XMLUni::fgXercesSchema, true);
    p_conf->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
    // Code taken from xsd/cxx/xml/dom/parsing-source.txx
    if (!rProps.schema_location().empty())
    {
        xml::string locn(rProps.schema_location());
        const void* p_locn(locn.c_str());
        p_conf->setParameter(XMLUni::fgXercesSchemaExternalSchemaLocation,
                             const_cast<void*>(p_locn));
    }
    if (!rProps.no_namespace_schema_location().empty())
    {
        xml::string locn(rProps.no_namespace_schema_location());
        const void* p_locn(locn.c_str());

        p_conf->setParameter(XMLUni::fgXercesSchemaExternalNoNameSpaceSchemaLocation,
                             const_cast<void*>(p_locn));
    }

    // We will release the DOM document ourselves.
    p_conf->setParameter(XMLUni::fgXercesUserAdoptsDOMDocument, true);

    // Set error handler.
    xml::dom::bits::error_handler_proxy<char> ehp(rErrorHandler);
    p_conf->setParameter(XMLUni::fgDOMErrorHandler, &ehp);

#else // _XERCES_VERSION >= 30000
    // Same as above but for Xerces-C++ 2 series.
    xml::dom::auto_ptr<DOMBuilder> p_parser(p_impl->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0));

    p_parser->setFeature(XMLUni::fgDOMComments, false);
    p_parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
    p_parser->setFeature(XMLUni::fgDOMEntities, false);
    p_parser->setFeature(XMLUni::fgDOMNamespaces, true);
    p_parser->setFeature(XMLUni::fgDOMWhitespaceInElementContent, false);
    p_parser->setFeature(XMLUni::fgDOMValidation, true);
    p_parser->setFeature(XMLUni::fgXercesSchema, true);
    p_parser->setFeature(XMLUni::fgXercesSchemaFullChecking, false);
    p_parser->setFeature(XMLUni::fgXercesUserAdoptsDOMDocument, true);

    // Code taken from xsd/cxx/xml/dom/parsing-source.txx
    if (!rProps.schema_location().empty())
    {
        xml::string locn(rProps.schema_location());
        const void* p_locn(locn.c_str());
        p_parser->setProperty(XMLUni::fgXercesSchemaExternalSchemaLocation,
                              const_cast<void*>(p_locn));
    }

    if (!rProps.no_namespace_schema_location().empty())
    {
        xml::string locn(rProps.no_namespace_schema_location());
        const void* p_locn(locn.c_str());

        p_parser->setProperty(XMLUni::fgXercesSchemaExternalNoNameSpaceSchemaLocation,
                              const_cast<void*>(p_locn));
    }

    xml::dom::bits::error_handler_proxy<char> ehp(rErrorHandler);
    p_parser->setErrorHandler(&ehp);

#endif // _XERCES_VERSION >= 30000

    // Do the parse
    xml::dom::auto_ptr<DOMDocument> p_doc(p_parser->parseURI(rFileName.c_str()));

    if (ehp.failed())
    {
        p_doc.reset();
    }

    return p_doc;
}

#define COVERAGE_IGNORE
void XmlTools::PrintNode(const std::string& rMsg, xercesc::DOMNode* pNode, bool showChildren)
{
    std::string prefix = xsd::cxx::xml::transcode<char>(pNode->getPrefix());
    std::string name = xsd::cxx::xml::transcode<char>(pNode->getLocalName());
    std::string nsuri = xsd::cxx::xml::transcode<char>(pNode->getNamespaceURI());
    std::cout << rMsg << " " << pNode << " " << prefix << ":" << name << " in " << nsuri << std::endl;
    if (showChildren)
    {
        for (xercesc::DOMNode* p_node = pNode->getFirstChild();
             p_node != NULL;
             p_node = p_node->getNextSibling())
        {
            std::cout << "     child type " << p_node->getNodeType();
            PrintNode("", p_node, false);
        }
        xercesc::DOMNamedNodeMap* p_attrs = pNode->getAttributes();
        if (p_attrs)
        {
            for (XMLSize_t i=0; i<p_attrs->getLength(); i++)
            {
                 xercesc::DOMNode* p_attr = p_attrs->item(i);
                 std::string value = xsd::cxx::xml::transcode<char>(p_attr->getNodeValue());
                 PrintNode("     attr (" + value + ")", p_attr, false);
            }
        }
    }
}
#undef COVERAGE_IGNORE

xercesc::DOMElement* XmlTools::SetNamespace(xercesc::DOMDocument* pDocument,
                                            xercesc::DOMElement* pElement,
                                            const XMLCh* pNamespace)
{
    using namespace xercesc;

    //PrintNode("Renaming", pElement, true);
    DOMNamedNodeMap* p_orig_attrs = pElement->getAttributes();
    std::vector<std::string> attr_values;
    if (p_orig_attrs)
    {
        for (XMLSize_t i=0; i<p_orig_attrs->getLength(); i++)
        {
            DOMNode* p_attr = p_orig_attrs->item(i);
            attr_values.push_back(xsd::cxx::xml::transcode<char>(p_attr->getNodeValue()));
        }
    }
    DOMElement* p_new_elt = static_cast<DOMElement*>(
        pDocument->renameNode(pElement, pNamespace, pElement->getLocalName()));
    //PrintNode("   to", p_new_elt, true);
    // Fix attributes - some get broken by the rename!
    if (p_orig_attrs)
    {
        DOMNamedNodeMap* p_new_attrs = p_new_elt->getAttributes();
        assert(p_new_attrs);
        assert(p_new_attrs == p_orig_attrs);
        assert(p_new_attrs->getLength() == attr_values.size());
        for (XMLSize_t i=0; i<p_new_attrs->getLength(); i++)
        {
            DOMNode* p_attr = p_new_attrs->item(i);
            p_attr->setNodeValue(X(attr_values[i]));
        }
    }
    //PrintNode("   after attr fix", p_new_elt, true);

    for (DOMNode* p_node = p_new_elt->getFirstChild();
         p_node != NULL;
         p_node = p_node->getNextSibling())
    {
        if (p_node->getNodeType() == DOMNode::ELEMENT_NODE)
        {
            p_node = SetNamespace(pDocument, static_cast<DOMElement*>(p_node), pNamespace);
        }
    }

    return p_new_elt;
}

xercesc::DOMElement* XmlTools::SetNamespace(xercesc::DOMDocument* pDocument,
                                            xercesc::DOMElement* pElement,
                                            const std::string& rNamespace)
{
    return SetNamespace(pDocument, pElement, X(rNamespace));
}

void XmlTools::FindElements(xercesc::DOMElement* pContextElement,
                            const std::vector<std::string>& rNames,
                            std::vector<xercesc::DOMElement*>& rResults,
                            unsigned depth)
{
    xercesc::DOMNodeList* p_child_elts = pContextElement->getElementsByTagName(X(rNames[depth]));
    unsigned num_children = p_child_elts->getLength();
    for (unsigned i=0; i<num_children; i++)
    {
        xercesc::DOMElement* p_child_elt = static_cast<xercesc::DOMElement*>(p_child_elts->item(i));
        if (depth == rNames.size() - 1)
        {
            rResults.push_back(p_child_elt);
        }
        else
        {
            FindElements(p_child_elt, rNames, rResults, depth+1);
        }
    }
}

std::vector<xercesc::DOMElement*> XmlTools::FindElements(xercesc::DOMElement* pContextElement,
                                                         const std::string& rPath)
{
    std::vector<xercesc::DOMElement*> results;
    std::vector<std::string> path;
    size_t start_pos = 0;
    size_t slash_pos = 0;
    while (slash_pos != std::string::npos)
    {
        slash_pos = rPath.find('/', start_pos);
        if (slash_pos == std::string::npos)
        {
            path.push_back(rPath.substr(start_pos));
        }
        else
        {
            path.push_back(rPath.substr(start_pos, slash_pos-start_pos));
        }
        start_pos = slash_pos + 1;
    }
    FindElements(pContextElement, path, results);
    return results;
}

void XmlTools::WrapContentInElement(xercesc::DOMDocument* pDocument,
                                    xercesc::DOMElement* pElement,
                                    const XMLCh* pNewElementLocalName)
{
    const XMLCh* p_namespace_uri = pElement->getNamespaceURI();
    const XMLCh* p_prefix = pElement->getPrefix();
    const XMLCh* p_qualified_name;
    if (p_prefix)
    {
#define COVERAGE_IGNORE
        // We can't actually cover this code, since versions of the parameters file which need this
        // transform didn't use a namespace, so can't have a namespace prefix!
        xercesc::QName qname(p_prefix, pNewElementLocalName, 0);
        p_qualified_name = qname.getRawName();
#undef COVERAGE_IGNORE
    }
    else
    {
        p_qualified_name = pNewElementLocalName;
    }
    xercesc::DOMElement* p_wrapper_elt = pDocument->createElementNS(p_namespace_uri, p_qualified_name);
    // Move all child nodes of pElement to be children of p_wrapper_elt
    xercesc::DOMNodeList* p_children = pElement->getChildNodes();
    for (unsigned i=0; i<p_children->getLength(); i++)
    {
        xercesc::DOMNode* p_child = pElement->removeChild(p_children->item(i));
        p_wrapper_elt->appendChild(p_child);
    }
    // Add the wrapper as the sole child of pElement
    pElement->appendChild(p_wrapper_elt);
}


std::string XmlTools::EscapeSpaces(const std::string& rPath)
{
    std::string escaped_path;
    for (std::string::const_iterator it = rPath.begin(); it != rPath.end(); ++it)
    {
        if (*it == ' ')
        {
            escaped_path += "%20";
        }
        else
        {
            escaped_path += *it;
        }
    }
    return escaped_path;
}
