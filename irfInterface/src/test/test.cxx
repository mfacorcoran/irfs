/**
 * @file test.cxx
 * @brief Test program for irfInterface.
 * @author J. Chiang
 *
 * $Header$
 */

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <algorithm>
#include <stdexcept>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "irfInterface/IrfsFactory.h"

#include "Aeff.h"
#include "Psf.h"
#include "Edisp.h"

using namespace irfInterface;

class irfInterfaceTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(irfInterfaceTests);

   CPPUNIT_TEST(test_addIrfs);
   CPPUNIT_TEST(test_create);
   CPPUNIT_TEST_EXCEPTION(test_creation_failure, std::invalid_argument);
   CPPUNIT_TEST(test_getIrfsNames);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void test_addIrfs();
   void test_create();
   void test_creation_failure();
   void test_getIrfsNames();

private:

   std::map<std::string, Irfs *> m_irfs;
   std::vector<std::string> m_irfNames;
   std::vector<std::string>::iterator m_name;
   IrfsFactory * m_irfsFactory;

};

void irfInterfaceTests::setUp() {
   m_irfNames.push_back("Moe");
   m_irfNames.push_back("Larry");
   m_irfNames.push_back("Curly");
   for (unsigned int i = 0; i < m_irfNames.size(); i++) {
      m_irfs[m_irfNames[i]] = new Irfs(new Aeff(), new Psf(), new Edisp(), i);
   }
   m_irfsFactory = IrfsFactory::instance();
}

void irfInterfaceTests::tearDown() {
   IrfsFactory::delete_instance();
   for (unsigned int i = 0; i < m_irfNames.size(); i++) {
      delete m_irfs[m_irfNames[i]];
   }
   m_irfNames.clear();
}

void irfInterfaceTests::test_addIrfs() {
   for (m_name = m_irfNames.begin(); m_name != m_irfNames.end(); ++m_name) {
      m_irfsFactory->addIrfs(*m_name, m_irfs[*m_name]->clone());
   }
}

void irfInterfaceTests::test_create() {
   test_addIrfs();
   for (unsigned int i = 0; i < m_irfNames.size(); i++) {
      Irfs * my_irfs = m_irfsFactory->create(m_irfNames[i]);
      CPPUNIT_ASSERT(my_irfs->irfID() == static_cast<int>(i));
   }
}

void irfInterfaceTests::test_creation_failure() {
   test_addIrfs();
   Irfs * my_irfs = m_irfsFactory->create("Shemp");
   (void)(my_irfs);
}

void irfInterfaceTests::test_getIrfsNames() {
   test_addIrfs();
   std::vector<std::string> names;
   m_irfsFactory->getIrfsNames(names);
   for (unsigned int i = 0; i < names.size(); i++) {
      m_name = std::find(m_irfNames.begin(), m_irfNames.end(), names[i]);
      CPPUNIT_ASSERT(m_name != m_irfNames.end());
   }
}

   
int main() {
   CppUnit::TextTestRunner runner;
   
   runner.addTest(irfInterfaceTests::suite());
    
   bool result = runner.run();
   if (result) {
      return 0;
   } else {
      return 1;
   }
}
