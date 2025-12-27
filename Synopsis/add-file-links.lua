function Pandoc(doc)
    local refs = doc.meta.references
    if refs then
        for i, ref in ipairs(refs) do
            if ref.file then
                -- This assumes your JabRef file field contains the path
                -- We wrap the title or a new [PDF] tag in a link
                local file_path = tostring(ref.file):gsub("^:", ""):gsub(":PDF$", "")
                ref.note = "[<a href='" .. file_path .. "'>Local PDF</a>]"
            end
        end
    end
    return doc
end