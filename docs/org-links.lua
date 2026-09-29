-- Turn relative links to other Org pages into Sphinx :doc: references.
--
-- The Org sources link to each other as [[file:page.org][text]] so they read
-- correctly on a forge. Pandoc would keep the .org target, which Sphinx does
-- not build, so a link to another page becomes a :doc: role on the same
-- relative path without the extension.

function Link(el)
  local target = el.target
  if target:match("^%a[%w+.-]*:") then
    return nil
  end
  local page = target:match("^(.-)%.org$")
  if not page then
    return nil
  end
  local text = pandoc.utils.stringify(el.content)
  return pandoc.RawInline("rst", ":doc:`" .. text .. " <" .. page .. ">`")
end
