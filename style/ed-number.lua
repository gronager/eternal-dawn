-- Eternal Dawn EPUB numbering, matching the PDF:
--   Parts (h1)            -> "Part I/II: <title>" (no section number)
--   Chapters (h2)         -> 0 (overture), 1..10, then "A" for the appendix
--   Sections (h3)         -> "<chap>.<n>"  (skipped if the chapter is the overture/unnumbered)
--   Deeper (h4+)          -> untouched (appendix Q-entries, paragraphs)
--   References/Bibliography & any \chapter*/\section* (unnumbered) -> left alone
local part, chap, sec, in_appendix, cur = 0, -1, 0, false, "0"
local roman = {"I","II","III","IV","V","VI"}

local function unnumbered(el)
  for _, c in ipairs(el.classes) do if c == "unnumbered" then return true end end
  return false
end
local function prefix(el, label)
  table.insert(el.content, 1, pandoc.Space())
  table.insert(el.content, 1, pandoc.Str(label))
end

function Header(el)
  local title = pandoc.utils.stringify(el)
  if title:find("References") or title:find("Bibliography") then return el end
  if el.level == 1 then
    part = part + 1
    table.insert(el.content, 1, pandoc.Str("Part " .. (roman[part] or part) .. ": "))
    return el
  elseif el.level == 2 then
    if unnumbered(el) then return el end
    sec = 0
    if in_appendix or title:find("Open Questions") then
      in_appendix = true; cur = "A"
    else
      chap = chap + 1; cur = tostring(chap)
    end
    prefix(el, cur)
  elseif el.level == 3 then
    if unnumbered(el) then return el end
    sec = sec + 1
    prefix(el, cur .. "." .. sec)
  end
  return el
end
