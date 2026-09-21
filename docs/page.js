function showPage(pageId) {
   document.querySelectorAll(".page").forEach(page => {
      page.classList.remove('active');
   });
   
   const target = document.getElementById(pageId);
   
   if (!target) {
      console.error("Page not found:", pageId);
      return;
   }
   
   target.classList.add("active");
   
   window.scrollTo({
      top: 0,
      behavior: 'auto'
   });
}

function showSection(pageId, sectionId) {
   showPage(pageId);
   const section = document.getElementById(sectionId);
   
   if (!section) {
      console.error("Section not found:", sectionId);
      return;
   }
   
   setTimeout(() => {
      section.scrollIntoView({
         behavior: "smooth",
         block: "start"
      });
   }, 0);
}